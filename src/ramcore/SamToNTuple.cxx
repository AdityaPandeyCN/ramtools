#include "ramcore/SamToNTuple.h"
#include "ramcore/SamParser.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <ROOT/RNTupleWriteOptions.hxx>
#include <TStopwatch.h>
#include <TList.h>
#include <TNamed.h>
#include <TFile.h>

#include <ROOT/RNTupleFillContext.hxx>
#include <ROOT/RNTupleParallelWriter.hxx>
#include <TROOT.h>

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <deque>
#include <exception>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace {

// Adds one header line to the list and registers the reference of an @SQ line,
// so that reference ids follow header order.
void HandleHeaderLine(TList &headers, const std::string &tag, const std::string &content)
{
   headers.Add(new TNamed(tag.c_str(), content.c_str()));

   if (tag == "@SQ") {
      size_t sn_pos = content.find("SN:");
      if (sn_pos != std::string::npos) {
         sn_pos += 3;
         size_t tab_pos = content.find('\t', sn_pos);
         std::string ref_name =
            content.substr(sn_pos, tab_pos != std::string::npos ? tab_pos - sn_pos : std::string::npos);
         RAMNTupleRecord::GetRnameRefs()->GetRefId(ref_name);
      }
   }
}

// Copies every field of a parsed SAM record into the RAM record except the two
// reference ids, which the caller resolves (the parallel converter caches them
// per thread).
void FillRecordFields(const ramcore::SamRecord &sam_record, RAMNTupleRecord &rec, uint32_t quality_policy)
{
   rec.SetBit(quality_policy);
   rec.SetQNAME(sam_record.qname);
   rec.SetFLAG(sam_record.flag);
   rec.SetPOS(sam_record.pos);
   rec.SetMAPQ(sam_record.mapq);
   rec.SetCIGAR(sam_record.cigar);
   rec.SetPNEXT(sam_record.pnext);
   rec.SetTLEN(sam_record.tlen);
   rec.SetSEQ(sam_record.seq);
   rec.SetQUAL(sam_record.qual);

   rec.ResetNOPT();
   for (const auto &opt : sam_record.optional_fields)
      rec.SetOPT(opt);
}

} // namespace

void samtoramntuple(const char *datafile, const char *treefile, bool split, bool cache, int compression_algorithm,
                    uint32_t quality_policy)
{
    TStopwatch stopwatch;
    stopwatch.Start();

    auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
    if (!rootFile || !rootFile->IsOpen()) {
        printf("Failed to create RAM file %s\n", treefile);
        return;
    }

    RAMNTupleRecord::InitializeRefs();

    auto model = RAMNTupleRecord::MakeModel();

    ROOT::RNTupleWriteOptions writeOptions;
    writeOptions.SetCompression(compression_algorithm);
    writeOptions.SetMaxUnzippedPageSize(64000);

    auto writer = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *rootFile, writeOptions);
    auto defaultEntry = writer->GetModel().CreateEntry();
    auto recordPtr = defaultEntry->GetPtr<RAMNTupleRecord>("record");

    TList headers;
    headers.SetName("headers");

    ramcore::SamParser parser;

    auto header_callback = [&headers](const std::string &tag, const std::string &content) {
       HandleHeaderLine(headers, tag, content);
    };

    auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t) {
       FillRecordFields(sam_record, *recordPtr, quality_policy);
       recordPtr->SetREFID(sam_record.rname);
       recordPtr->SetREFNEXT(sam_record.rnext);

       RAMNTupleRecord::NoteRefSpan(recordPtr->GetRefSpan());
       RAMNTupleRecord::NotePlacement(recordPtr->GetREFID(), recordPtr->GetPOS() - 1);
       writer->Fill(*defaultEntry);
    };

    if (!parser.ParseFile(datafile, header_callback, record_callback)) {
        printf("Failed to parse SAM file %s\n", datafile);
        return;
    }

    writer.reset();

    // Region queries can only seek on a sorted file; the file records which it is.
    if (!RAMNTupleRecord::IsCoordinateSorted())
       fprintf(stderr, "%s is not in coordinate order; region queries will read it in full.\n", datafile);
    RAMNTupleRecord::WriteAllRefs(*rootFile);

    // One key for the list; without kSingleKey every line is written as its own
    // key and no reader can get the header back in order.
    headers.Write("headers", TObject::kSingleKey);
    rootFile->Close();

    printf("\nRAM file created: %s\n", treefile);
    printf("Number of entries: %zu\n", parser.GetRecordsProcessed());

    RAMNTupleRecord::GetRnameRefs()->Print();
    RAMNTupleRecord::GetRnextRefs()->Print();

    printf("\nProcessed %zu SAM headers\n", parser.GetLinesProcessed() - parser.GetRecordsProcessed());
    printf("Processed %zu SAM records\n\n", parser.GetRecordsProcessed());

    stopwatch.Print();
}

namespace {

// One output file, held open while the input streams past it.
struct ChromosomeWriter {
   std::unique_ptr<TFile> file{};
   std::unique_ptr<ROOT::RNTupleWriter> writer{};
   std::unique_ptr<ROOT::REntry> entry{};
   std::shared_ptr<RAMNTupleRecord> record{};
   int64_t rows = 0;
   int32_t last_pos = -1;
   bool sorted = true;
};

} // namespace

void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy)
{
   RAMNTupleRecord::InitializeRefs();

   std::map<std::string, ChromosomeWriter> writers;
   TList headers;
   headers.SetName("headers");
   headers.SetOwner(true);

   auto header_callback = [&](const std::string &tag, const std::string &content) {
      HandleHeaderLine(headers, tag, content);
   };

   auto open_writer = [&](const std::string &chr) -> ChromosomeWriter & {
      auto it = writers.find(chr);
      if (it != writers.end())
         return it->second;

      ChromosomeWriter &cw = writers[chr];
      std::string filename{output_prefix};
      filename += "_";
      filename += chr;
      filename += ".root";
      cw.file.reset(TFile::Open(filename.c_str(), "RECREATE"));
      if (!cw.file || !cw.file->IsOpen())
         throw std::runtime_error("cannot create " + filename);

      ROOT::RNTupleWriteOptions writeOptions;
      writeOptions.SetCompression(compression_algorithm);
      writeOptions.SetMaxUnzippedPageSize(64000);
      // Every chromosome's buffers are open at once, so keep the clusters small.
      writeOptions.SetApproxZippedClusterSize(8 * 1024 * 1024);

      cw.writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *cw.file, writeOptions);
      cw.entry = cw.writer->GetModel().CreateEntry();
      cw.record = cw.entry->GetPtr<RAMNTupleRecord>("record");
      return cw;
   };

   auto record_callback = [&](const ramcore::SamRecord &sam_record, size_t) {
      // A record with no reference has no chromosome file to go to.
      if (sam_record.rname == "*")
         return;

      ChromosomeWriter &cw = open_writer(sam_record.rname);
      RAMNTupleRecord &rec = *cw.record;

      FillRecordFields(sam_record, rec, quality_policy);
      rec.SetREFID(sam_record.rname);
      rec.SetREFNEXT(sam_record.rnext);

      RAMNTupleRecord::NoteRefSpan(rec.GetRefSpan());
      cw.writer->Fill(*cw.entry);
      cw.rows++;

      // One reference per file, so the order check is on the position alone.
      const int32_t pos = rec.GetPOS() - 1;
      if (pos < cw.last_pos)
         cw.sorted = false;
      cw.last_pos = pos;
   };

   ramcore::SamParser parser;
   if (!parser.ParseFile(datafile, header_callback, record_callback)) {
      printf("Failed to parse SAM file %s\n", datafile);
      return;
   }

   // The reference table and the longest span are only complete once the whole
   // input has been read, so every file is finished here.
   for (auto &[chr, cw] : writers) {
      cw.writer.reset();

      RAMNTupleRecord::SetCoordinateSorted(cw.sorted);
      if (!cw.sorted)
         fprintf(stderr, "%s: %s is not in coordinate order; region queries will read it in full.\n", datafile,
                 chr.c_str());
      RAMNTupleRecord::WriteAllRefs(*cw.file);

      cw.file->cd();
      headers.Write("headers", TObject::kSingleKey);
      cw.file->Close();

      printf("%s_%s.root: %lld records\n", output_prefix, chr.c_str(), static_cast<long long>(cw.rows));
   }
}

// ---------------------------------------------------------------------------
// Parallel conversion
//
// The main thread cuts the input into blocks of whole lines and hands them to
// worker threads through a bounded queue. Every worker owns an RNTupleFillContext
// of a shared RNTupleParallelWriter, so it parses, encodes and compresses its
// blocks without touching the other workers. The fill contexts stage their
// clusters instead of committing them; a worker commits the staged clusters of
// block n only once block n-1 has been committed, so the clusters, and with them
// the records, come out in input order.

namespace {

// Per-block summary of what the order check and the longest span need.
struct BlockOrder {
   bool has_placed = false;
   int32_t first_refid = -1;
   int32_t first_pos = -1;
   int32_t last_refid = -1;
   int32_t last_pos = -1;
   bool seen_unplaced = false;
   /// Placed records in non-decreasing order and none of them after an unplaced one.
   bool sorted = true;
   uint32_t max_span = 0;

   void Note(int32_t refid, int32_t pos)
   {
      if (refid < 0) {
         seen_unplaced = true;
         return;
      }
      if (seen_unplaced || (has_placed && (refid < last_refid || (refid == last_refid && pos < last_pos))))
         sorted = false;
      if (!has_placed) {
         has_placed = true;
         first_refid = refid;
         first_pos = pos;
      }
      last_refid = refid;
      last_pos = pos;
   }
};

// The same check as RAMNTupleRecord::NotePlacement, applied block by block in
// input order.
struct FileOrder {
   bool sorted = true;
   bool seen_unplaced = false;
   bool has_placed = false;
   int32_t last_refid = -1;
   int32_t last_pos = -1;
   uint32_t max_span = 0;

   void Add(const BlockOrder &b)
   {
      if (!b.sorted)
         sorted = false;
      if (b.has_placed) {
         if (seen_unplaced || (has_placed && (b.first_refid < last_refid ||
                                              (b.first_refid == last_refid && b.first_pos < last_pos))))
            sorted = false;
         has_placed = true;
         last_refid = b.last_refid;
         last_pos = b.last_pos;
      }
      seen_unplaced = seen_unplaced || b.seen_unplaced;
      max_span = std::max(max_span, b.max_span);
   }
};

struct Block {
   size_t seq = 0;        ///< Position in the input, from 0.
   size_t first_line = 0; ///< 1-based number of the block's first line, for warnings.
   std::vector<char> data;
};

// Reads the input in blocks of whole lines. A line longer than a block is read
// in full by growing the block.
class BlockReader {
   FILE *fFile;
   size_t fBlockBytes;
   std::vector<char> fCarry;
   bool fEof = false;

public:
   BlockReader(FILE *file, size_t block_bytes) : fFile(file), fBlockBytes(std::max<size_t>(block_bytes, 1)) {}

   /// Returns false once the input is exhausted. Every line in \p out ends in
   /// '\n' except possibly the last line of the file.
   bool Next(std::vector<char> &out)
   {
      out.swap(fCarry);
      fCarry.clear();
      while (!fEof) {
         const size_t old = out.size();
         out.resize(old + fBlockBytes);
         const size_t n = fread(out.data() + old, 1, fBlockBytes, fFile);
         out.resize(old + n);
         if (n < fBlockBytes) {
            if (ferror(fFile))
               throw std::runtime_error("read error on the SAM input");
            fEof = true;
         }
         const auto *end = out.data() + out.size();
         const auto *nl = static_cast<const char *>(memrchr(out.data(), '\n', out.size()));
         if (nl) {
            fCarry.assign(nl + 1, end);
            out.resize(static_cast<size_t>(nl + 1 - out.data()));
            return true;
         }
         // No newline yet: the line is longer than a block, keep reading.
      }
      return !out.empty();
   }
};

// Blocks travel from the reader to the workers through here. The capacity
// bounds how far the reader runs ahead, and so the memory in flight.
class BlockQueue {
   std::mutex fMutex;
   std::condition_variable fNotEmpty;
   std::condition_variable fNotFull;
   std::deque<Block> fBlocks;
   size_t fCapacity;
   bool fClosed = false;

public:
   explicit BlockQueue(size_t capacity) : fCapacity(std::max<size_t>(capacity, 1)) {}

   void Push(Block &&block)
   {
      std::unique_lock<std::mutex> lock(fMutex);
      fNotFull.wait(lock, [&] { return fBlocks.size() < fCapacity || fClosed; });
      if (fClosed)
         return;
      fBlocks.push_back(std::move(block));
      fNotEmpty.notify_one();
   }

   /// Returns false when the queue is closed and drained.
   bool Pop(Block &block)
   {
      std::unique_lock<std::mutex> lock(fMutex);
      fNotEmpty.wait(lock, [&] { return !fBlocks.empty() || fClosed; });
      if (fBlocks.empty())
         return false;
      block = std::move(fBlocks.front());
      fBlocks.pop_front();
      fNotFull.notify_one();
      return true;
   }

   /// Wakes everyone; Push() drops its block and Pop() returns false once drained.
   void Close()
   {
      const std::lock_guard<std::mutex> lock(fMutex);
      fClosed = true;
      fNotEmpty.notify_all();
      fNotFull.notify_all();
   }

   /// Close() and drop what is queued, for shutting down after an error.
   void Abort()
   {
      const std::lock_guard<std::mutex> lock(fMutex);
      fClosed = true;
      fBlocks.clear();
      fNotEmpty.notify_all();
      fNotFull.notify_all();
   }
};

// Remembers the last name looked up, so a sorted file resolves its reference
// almost always without taking the table's lock.
struct RefCache {
   std::string name;
   int id = -1;
   bool valid = false;

   int Lookup(RAMNTupleRefs &refs, const std::string &rname)
   {
      if (valid && rname == name)
         return id;
      id = refs.GetRefId(rname);
      name = rname;
      valid = true;
      return id;
   }
};

// What the workers report back, merged on the main thread in block order.
struct Progress {
   std::mutex mutex;
   std::condition_variable next_turn;
   size_t next_seq = 0; ///< The block whose clusters may be committed next.
   bool failed = false;
   std::exception_ptr error;
   FileOrder order;
   size_t records = 0;
   /// Header lines found among the records, in input order; the sequential
   /// converter accepts them anywhere, so this one does too.
   std::vector<std::pair<std::string, std::string>> late_headers;

   void Fail(std::exception_ptr e)
   {
      const std::lock_guard<std::mutex> lock(mutex);
      if (!failed) {
         failed = true;
         error = std::move(e);
      }
      next_turn.notify_all();
   }
};

// Encodes the records of one block into the worker's fill context and returns
// the block's order summary.
BlockOrder ProcessBlock(Block &block, ROOT::RNTupleFillContext &ctx, ROOT::REntry &entry, RAMNTupleRecord &rec,
                        uint32_t quality_policy, RefCache &rname_cache, RefCache &rnext_cache,
                        ramcore::SamRecord &sam_record, size_t &records,
                        std::vector<std::pair<std::string, std::string>> &late_headers)
{
   BlockOrder order;
   RAMNTupleRefs &rname_refs = *RAMNTupleRecord::GetRnameRefs();
   RAMNTupleRefs &rnext_refs = *RAMNTupleRecord::GetRnextRefs();

   char *cursor = block.data.data();
   char *const end = cursor + block.data.size();
   size_t line_number = block.first_line;
   while (cursor < end) {
      char *nl = static_cast<char *>(memchr(cursor, '\n', static_cast<size_t>(end - cursor)));
      char *line_end = nl ? nl : end;
      char *line = cursor;
      cursor = nl ? nl + 1 : end;

      while (line_end > line && (line_end[-1] == '\n' || line_end[-1] == '\r'))
         --line_end;
      *line_end = '\0';
      const size_t this_line = line_number++;

      if (line == line_end)
         continue;

      if (line[0] == '@') {
         char *tab = strchr(line, '\t');
         if (tab) {
            *tab = '\0';
            late_headers.emplace_back(line, tab + 1);
         } else {
            late_headers.emplace_back(line, "");
         }
         continue;
      }

      sam_record.Clear();
      if (!ramcore::SamParser::ParseRecord(line, sam_record, this_line))
         continue;

      FillRecordFields(sam_record, rec, quality_policy);
      rec.refid = rname_cache.Lookup(rname_refs, sam_record.rname);
      rec.refnext = rnext_cache.Lookup(rnext_refs, sam_record.rnext);

      order.max_span = std::max(order.max_span, rec.GetRefSpan());
      order.Note(rec.refid, rec.pos);
      ctx.Fill(entry);
      records++;
   }
   return order;
}

void WorkerMain(BlockQueue &queue, Progress &progress, std::shared_ptr<ROOT::RNTupleFillContext> ctx,
                uint32_t quality_policy)
{
   try {
      auto entry = ctx->CreateEntry();
      auto rec = entry->GetPtr<RAMNTupleRecord>("record");
      RefCache rname_cache;
      RefCache rnext_cache;
      ramcore::SamRecord sam_record;

      Block block;
      while (queue.Pop(block)) {
         size_t records = 0;
         std::vector<std::pair<std::string, std::string>> late_headers;
         const BlockOrder order = ProcessBlock(block, *ctx, *entry, *rec, quality_policy, rname_cache, rnext_cache,
                                               sam_record, records, late_headers);
         // Writes the block's pages to the file and stages the cluster; the
         // logical append below is what has to wait for its turn.
         ctx->FlushCluster();

         std::unique_lock<std::mutex> lock(progress.mutex);
         progress.next_turn.wait(lock, [&] { return progress.next_seq == block.seq || progress.failed; });
         if (progress.failed)
            return;
         ctx->CommitStagedClusters();
         progress.order.Add(order);
         progress.records += records;
         progress.late_headers.insert(progress.late_headers.end(), late_headers.begin(), late_headers.end());
         progress.next_seq++;
         progress.next_turn.notify_all();
      }
   } catch (...) {
      progress.Fail(std::current_exception());
      queue.Abort();
   }
}

// Consumes the header lines at the front of \p data, returning the offset of
// the first record line or npos when the whole block is header.
size_t ConsumeHeader(std::vector<char> &data, TList &headers, size_t &lines)
{
   size_t offset = 0;
   while (offset < data.size()) {
      const auto *nl = static_cast<const char *>(memchr(data.data() + offset, '\n', data.size() - offset));
      const size_t next = nl ? static_cast<size_t>(nl - data.data()) + 1 : data.size();
      size_t line_end = nl ? static_cast<size_t>(nl - data.data()) : data.size();
      while (line_end > offset && (data[line_end - 1] == '\r' || data[line_end - 1] == '\n'))
         --line_end;
      if (line_end == offset) {
         // Empty line: skipped, as the sequential parser skips it.
         offset = next;
         lines++;
         continue;
      }
      if (data[offset] != '@')
         return offset;
      std::string line(data.data() + offset, line_end - offset);
      const size_t tab = line.find('\t');
      if (tab != std::string::npos)
         HandleHeaderLine(headers, line.substr(0, tab), line.substr(tab + 1));
      else
         HandleHeaderLine(headers, line, "");
      offset = next;
      lines++;
   }
   return std::string::npos;
}

} // namespace

void samtoramntuple_parallel(const char *datafile, const char *treefile, int compression_algorithm,
                             uint32_t quality_policy, int threads, size_t block_bytes)
{
   TStopwatch stopwatch;
   stopwatch.Start();

   threads = std::max(threads, 1);

   std::unique_ptr<FILE, int (*)(FILE *)> input(fopen(datafile, "r"), fclose);
   if (!input) {
      printf("Failed to parse SAM file %s\n", datafile);
      return;
   }

   auto rootFile = std::unique_ptr<TFile>(TFile::Open(treefile, "RECREATE"));
   if (!rootFile || !rootFile->IsOpen()) {
      printf("Failed to create RAM file %s\n", treefile);
      return;
   }

   // The workers construct records, look reference names up and stream through
   // ROOT's type system at the same time.
   ROOT::EnableThreadSafety();
   RAMNTupleRecord::InitializeRefs();

   TList headers;
   headers.SetName("headers");

   // The header comes first, so it is read here before any worker starts: the
   // @SQ lines define the reference ids.
   BlockReader reader(input.get(), block_bytes);
   size_t lines = 0;
   Block first;
   bool have_records = false;
   {
      std::vector<char> data;
      while (reader.Next(data)) {
         const size_t offset = ConsumeHeader(data, headers, lines);
         if (offset == std::string::npos) {
            data.clear();
            continue;
         }
         first.data.assign(data.begin() + static_cast<std::ptrdiff_t>(offset), data.end());
         have_records = true;
         break;
      }
   }

   // RNEXT ids in header order rather than in order of first appearance, which
   // would depend on which worker gets there first. "=" is the common case and
   // keeps id 0, as in the sequential converter.
   {
      RAMNTupleRefs &rnext = *RAMNTupleRecord::GetRnextRefs();
      rnext.GetRefId("=");
      for (const auto &name : RAMNTupleRecord::GetRnameRefs()->GetRefs())
         rnext.GetRefId(name);
   }

   auto model = ROOT::RNTupleModel::CreateBare();
   model->MakeField<RAMNTupleRecord>("record");

   ROOT::RNTupleWriteOptions writeOptions;
   writeOptions.SetCompression(compression_algorithm);
   writeOptions.SetMaxUnzippedPageSize(64000);
   // The parallel writer needs buffered writing (the default); said explicitly
   // because it is a requirement, not a tuning choice.
   writeOptions.SetUseBufferedWrite(true);

   Progress progress;
   {
      auto writer = ROOT::RNTupleParallelWriter::Append(std::move(model), "RAM", *rootFile, writeOptions);

      // One block per worker in flight plus one per worker queued bounds the
      // memory at about 2 * threads * block_bytes.
      BlockQueue queue(static_cast<size_t>(threads));
      std::vector<std::shared_ptr<ROOT::RNTupleFillContext>> contexts;
      std::vector<std::thread> workers;
      for (int i = 0; i < threads; i++) {
         auto ctx = writer->CreateFillContext();
         ctx->EnableStagedClusterCommitting();
         contexts.push_back(ctx);
         workers.emplace_back(WorkerMain, std::ref(queue), std::ref(progress), ctx, quality_policy);
      }

      if (have_records) {
         size_t seq = 0;
         Block block = std::move(first);
         do {
            block.seq = seq++;
            block.first_line = lines + 1;
            lines += static_cast<size_t>(std::count(block.data.begin(), block.data.end(), '\n'));
            if (!block.data.empty() && block.data.back() != '\n')
               lines++;
            queue.Push(std::move(block));
            block = Block{};
            {
               const std::lock_guard<std::mutex> lock(progress.mutex);
               if (progress.failed)
                  break;
            }
         } while (reader.Next(block.data));
      }
      queue.Close();
      for (auto &w : workers)
         w.join();

      // Every context has to be gone before the writer commits the dataset.
      contexts.clear();
      if (progress.failed) {
         writer.reset();
         rootFile->Close();
         std::rethrow_exception(progress.error);
      }
      writer.reset();
   }

   for (const auto &[tag, content] : progress.late_headers)
      HandleHeaderLine(headers, tag, content);

   RAMNTupleRecord::SetCoordinateSorted(progress.order.sorted);
   RAMNTupleRecord::NoteRefSpan(progress.order.max_span);
   if (!progress.order.sorted)
      fprintf(stderr, "%s is not in coordinate order; region queries will read it in full.\n", datafile);
   RAMNTupleRecord::WriteAllRefs(*rootFile);

   headers.Write("headers", TObject::kSingleKey);
   rootFile->Close();

   printf("\nRAM file created: %s\n", treefile);
   printf("Number of entries: %zu\n", progress.records);

   RAMNTupleRecord::GetRnameRefs()->Print();
   RAMNTupleRecord::GetRnextRefs()->Print();

   printf("\nProcessed %d SAM headers\n", headers.GetSize());
   printf("Processed %zu SAM records with %d threads\n\n", progress.records, threads);

   stopwatch.Print();
}
