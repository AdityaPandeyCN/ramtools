#include "ramcore/QualityBlocks.h"
#include "ramcore/SeqBlocks.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

const char *const kFile = "seq_blocks_test.root";
const std::string kBases = "ACGT";

std::string RandomBases(std::mt19937 &rng, std::size_t n)
{
   std::string s(n, 'A');
   for (auto &c : s)
      c = kBases[rng() % 4];
   return s;
}

struct Read {
   std::string rname; ///< "*" for an unplaced read
   int32_t pos;       ///< 1-based, as in SAM
   std::string cigar;
   std::string seq;
};

/// Overlapping reads from one random reference, with some sequencing errors.
std::vector<Read> PileupReads(std::size_t count)
{
   std::mt19937 rng(7);
   const std::string reference = RandomBases(rng, /*n=*/5000);
   std::vector<Read> reads;
   for (std::size_t i = 0; i < count; i++) {
      const std::size_t start = (i * 37) % (reference.size() - 100);
      std::string seq = reference.substr(start, 100);
      if (rng() % 4 == 0)
         seq[rng() % seq.size()] = kBases[rng() % 4];
      reads.push_back({"chr1", static_cast<int32_t>(start + 1), "100M", seq});
   }
   return reads;
}

void Write(const std::vector<Read> &reads, std::size_t block)
{
   RAMNTupleRecord::InitializeRefs();
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      QualityBlockWriter out(
         writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
         [&writer](ROOT::REntry &e) { writer->Fill(e); }, block);
      for (const auto &read : reads) {
         RAMNTupleRecord &rec = out.Record();
         rec.SetBit(RAMNTupleRecord::kPhred33);
         rec.SetQNAME("r");
         rec.SetFLAG(read.rname == "*" ? 4 : 0);
         rec.SetRNAME(read.rname);
         rec.SetPOS(read.pos);
         rec.SetCIGAR(read.cigar);
         rec.SetSEQ(read.seq);
         rec.SetQUAL("*");
         out.Add();
      }
      out.Finish();
      RAMNTupleRecord::SetQualBlockEnds(out.TakeBlockEnds());
   }
   RAMNTupleRecord::WriteAllRefs(*file);
   file->Close();
}

void ExpectSeqs(const std::vector<Read> &reads, const std::vector<std::size_t> &order)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   SeqBlockReader seqs(*reader);
   for (const std::size_t row : order)
      EXPECT_EQ(seqs.Get(view(row), row), reads[row].seq) << "row " << row;
}

std::vector<std::size_t> InOrder(std::size_t n)
{
   std::vector<std::size_t> order(n);
   for (std::size_t i = 0; i < n; i++)
      order[i] = i;
   return order;
}

class SeqBlocksTest : public ::testing::Test {
protected:
   void TearDown() override { std::remove(kFile); }
};

TEST(SeqBlockCodecTest, CodesOnlyACGTN)
{
   EXPECT_TRUE(SeqBlockCodec::Codable("ACGTN"));
   EXPECT_TRUE(SeqBlockCodec::Codable("NNNN"));
   EXPECT_FALSE(SeqBlockCodec::Codable(""));
   EXPECT_FALSE(SeqBlockCodec::Codable("*"));
   EXPECT_FALSE(SeqBlockCodec::Codable("ACGR"));
   EXPECT_FALSE(SeqBlockCodec::Codable("AC=T"));
   EXPECT_FALSE(SeqBlockCodec::Codable("acgt"));
}

TEST(SeqBlockCodecTest, EveryReadComesBackExactly)
{
   std::mt19937 rng(3);
   std::vector<SeqBlockRead> reads;
   std::vector<std::string> seqs;
   auto add = [&](int32_t refid, int32_t pos, std::vector<uint32_t> cigar, std::string seq) {
      reads.push_back({refid, pos, std::move(cigar)});
      seqs.push_back(std::move(seq));
   };
   const std::string ref = RandomBases(rng, /*n=*/400);
   const auto m = [](uint32_t n) { return n << 4; }; // BAM CIGAR op M
   add(0, 100, {m(50)}, ref.substr(0, 50));
   add(0, 110, {m(50)}, ref.substr(10, 50));
   add(0, 110, {m(50)}, ref.substr(10, 20) + "N" + ref.substr(31, 29)); // one N
   add(0, 120, {m(50)}, "NNNNN" + ref.substr(25, 40) + "NNNNN");        // N at both ends
   add(0, 120, {m(50)}, std::string(50, 'N'));                          // only N
   add(0, 130, {(3U << 4) | 4U, m(47)}, "GGG" + ref.substr(33, 47));    // soft clip
   add(0, 130, {m(10), (2U << 4) | 1U, m(38)}, ref.substr(30, 50));     // insertion
   add(0, 140, {m(20), (5U << 4) | 2U, m(30)}, ref.substr(40, 50));     // deletion
   add(-1, -1, {}, RandomBases(rng, /*n=*/77));                         // no CIGAR: length coded
   add(-1, -1, {}, RandomBases(rng, /*n=*/77));                         // same length as before
   add(0, 150, {m(50)}, ref.substr(50, 60));                            // longer than its CIGAR
   add(0, 150, {m(50)}, ref.substr(50, 3));                             // shorter than its CIGAR
   add(0, 150, {}, "A");                                                // one base
   add(1, 150, {m(50)}, ref.substr(0, 50));                             // other reference
   add(0, 0, {m(100000)}, RandomBases(rng, /*n=*/100000));              // long read
   for (int i = 0; i < 200; i++)                                        // repeats the match model follows
      add(-1, -1, {}, ref.substr(static_cast<std::size_t>(i), 100));

   const auto bytes = SeqBlockCodec::Encode(reads, seqs);
   EXPECT_EQ(SeqBlockCodec::Decode(bytes, reads), seqs);
}

TEST(SeqBlockCodecTest, RejectsBlocksOfAnotherFormat)
{
   const std::vector<SeqBlockRead> reads{{-1, -1, {}}};
   auto bytes = SeqBlockCodec::Encode(reads, {"ACGT"});
   bytes[0]++;
   EXPECT_THROW(SeqBlockCodec::Decode(bytes, reads), std::runtime_error);
   EXPECT_THROW(SeqBlockCodec::Encode(reads, {"ACGR"}), std::invalid_argument);
}

TEST_F(SeqBlocksTest, OverlappingReadsCodeInFarFewerBytes)
{
   const std::vector<Read> reads = PileupReads(/*count=*/2000);
   std::vector<SeqBlockRead> in;
   std::vector<std::string> seqs;
   RAMNTupleRecord::InitializeRefs();
   for (const auto &r : reads) {
      RAMNTupleRecord rec;
      rec.SetCIGAR(r.cigar);
      in.push_back({0, r.pos - 1, rec.cigar});
      seqs.push_back(r.seq);
   }
   const auto bytes = SeqBlockCodec::Encode(in, seqs);
   // 200,000 bases; the pileup predicts nearly all of them.
   EXPECT_LT(bytes.size(), 10000U);
   EXPECT_EQ(SeqBlockCodec::Decode(bytes, in), seqs);
}

TEST_F(SeqBlocksTest, ReadsEveryRowInAnyOrder)
{
   std::vector<Read> reads = PileupReads(/*count=*/1000);
   // Reads the model does not code stay in the record.
   reads[3].seq = "ACGRYACGT";
   reads[500].seq = "*";
   reads[501].rname = "*";
   reads[501].pos = 0;
   reads[501].cigar = "*";
   reads[999].seq = "NNNNACGTNN";
   Write(reads, /*block=*/64);

   std::vector<std::size_t> order = InOrder(reads.size());
   ExpectSeqs(reads, order);
   std::mt19937 rng(11);
   std::shuffle(order.begin(), order.end(), rng);
   ExpectSeqs(reads, order);

   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   EXPECT_FALSE(view(3).TestBit(RAMNTupleRecord::kSeqInBlock));
   EXPECT_EQ(view(3).seq, "ACGRYACGT");
   EXPECT_FALSE(view(500).TestBit(RAMNTupleRecord::kSeqInBlock));
   EXPECT_TRUE(view(501).TestBit(RAMNTupleRecord::kSeqInBlock));
   EXPECT_TRUE(view(999).TestBit(RAMNTupleRecord::kSeqInBlock));
   EXPECT_TRUE(view(999).seq.empty());
}

TEST_F(SeqBlocksTest, BlockWithNoCodableReadStoresNothing)
{
   const std::vector<Read> reads{{"chr1", 1, "4M", "ACGR"}, {"chr1", 2, "1M", "*"}};
   Write(reads, /*block=*/10);
   ExpectSeqs(reads, InOrder(reads.size()));
}

} // namespace
