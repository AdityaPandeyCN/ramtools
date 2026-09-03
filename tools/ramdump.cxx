//
// ramdump -- write a RAM file back out as SAM, the way `samtools view` does.
//
// Having a text dump is what makes the RAM format checkable: the output of
// `ramdump -h file.ram` should be identical to `samtools view -h` on the
// alignments the file was built from, and `ramdump -c region` should agree with
// `samtools view -c region`.
//

#include "ramcore/RAMNTupleView.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleReader.hxx>
#include <TFile.h>
#include <TList.h>
#include <TNamed.h>

#include <charconv>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>

namespace {

// SAM is line-oriented text over millions of records, so every field goes
// through one growing buffer that is handed to fwrite in big blocks. Formatting
// straight into std::ostream costs more than the query itself.
class SamWriter {
public:
   explicit SamWriter(FILE *out) : fOut(out) { fBuf.reserve(kFlushAt + 4096); }
   ~SamWriter() { Flush(); }

   SamWriter(const SamWriter &) = delete;
   SamWriter &operator=(const SamWriter &) = delete;

   void Str(const std::string &s) { fBuf.append(s); }
   void Str(const char *s) { fBuf.append(s); }
   void Char(char c) { fBuf.push_back(c); }

   /// Direct access for fields that decode straight into the buffer.
   std::string &Buffer() { return fBuf; }

   void Int(long long v)
   {
      char tmp[24];
      auto res = std::to_chars(tmp, tmp + sizeof(tmp), v);
      fBuf.append(tmp, res.ptr - tmp);
   }

   void EndLine()
   {
      fBuf.push_back('\n');
      if (fBuf.size() >= kFlushAt)
         Flush();
   }

   void Flush()
   {
      if (!fBuf.empty()) {
         fwrite(fBuf.data(), 1, fBuf.size(), fOut);
         fBuf.clear();
      }
   }

private:
   static constexpr size_t kFlushAt = 1u << 20;
   FILE *fOut;
   std::string fBuf;
};

// The converter stores each header line as a TNamed holding the tag ("@SQ") and
// everything after the first tab, so the original line is tag + '\t' + content.
// A tag-only line such as a bare "@CO" has empty content and gets no tab back.
void WriteHeader(const char *file, SamWriter &out)
{
   std::unique_ptr<TFile> f(TFile::Open(file, "READ"));
   if (!f || f->IsZombie())
      return;

   auto *headers = f->Get<TList>("headers");
   if (!headers)
      return;

   TIter next(headers);
   while (auto *obj = next()) {
      auto *named = dynamic_cast<TNamed *>(obj);
      if (!named)
         continue;
      out.Str(named->GetName());
      const char *content = named->GetTitle();
      if (content && *content) {
         out.Char('\t');
         out.Str(content);
      }
      out.EndLine();
   }
}

void WriteRecord(const RAMNTupleRecord &rec, SamWriter &out)
{
   out.Str(rec.GetQNAME());
   out.Char('\t');
   out.Int(rec.GetFLAG());
   out.Char('\t');
   out.Str(rec.GetRNAME());
   out.Char('\t');
   out.Int(rec.GetPOS());
   out.Char('\t');
   out.Int(rec.GetMAPQ());
   out.Char('\t');
   rec.AppendCIGAR(out.Buffer());
   out.Char('\t');
   out.Str(rec.GetRNEXT());
   out.Char('\t');
   out.Int(rec.GetPNEXT());
   out.Char('\t');
   out.Int(rec.GetTLEN());
   out.Char('\t');
   rec.AppendSEQ(out.Buffer());
   out.Char('\t');
   rec.AppendQUAL(out.Buffer());

   for (const auto &tag : rec.GetTags()) {
      out.Char('\t');
      out.Str(tag);
   }
   out.EndLine();
}

void Usage(const char *argv0)
{
   std::cerr << "Usage: " << argv0 << " [options] <in.ram> [region]\n"
             << "Options:\n"
             << "  -h        include the header in the output\n"
             << "  -H        print the header only\n"
             << "  -c        print only the number of matching records\n"
             << "  -f INT    only records with all of these FLAG bits set\n"
             << "  -F INT    only records with none of these FLAG bits set\n"
             << "  -o FILE   write to FILE instead of stdout\n"
             << "\nRegion is rname[:start[-end]], 1-based and inclusive, as in samtools.\n";
}

bool ParseFlag(const char *text, uint16_t &out)
{
   char *end = nullptr;
   const long value = std::strtol(text, &end, 0); // 0: accepts 0x900 as well as 2304
   if (end == text || *end != '\0' || value < 0 || value > 0xFFFF)
      return false;
   out = static_cast<uint16_t>(value);
   return true;
}

} // namespace

int main(int argc, char *argv[])
{
   bool withHeader = false;
   bool headerOnly = false;
   bool countOnly = false;
   uint16_t requireFlags = 0;
   uint16_t excludeFlags = 0;
   const char *outPath = nullptr;
   const char *file = nullptr;
   const char *region = "";

   for (int i = 1; i < argc; i++) {
      const std::string arg = argv[i];
      if (arg == "-h") {
         withHeader = true;
      } else if (arg == "-H") {
         headerOnly = true;
      } else if (arg == "-c") {
         countOnly = true;
      } else if (arg == "-f" || arg == "-F") {
         if (i + 1 >= argc) {
            std::cerr << "ramdump: " << arg << " needs a value\n";
            return 1;
         }
         uint16_t value = 0;
         if (!ParseFlag(argv[++i], value)) {
            std::cerr << "ramdump: invalid FLAG value '" << argv[i] << "'\n";
            return 1;
         }
         (arg == "-f" ? requireFlags : excludeFlags) = value;
      } else if (arg == "-o") {
         if (i + 1 >= argc) {
            std::cerr << "ramdump: -o needs a file name\n";
            return 1;
         }
         outPath = argv[++i];
      } else if (arg == "--help") {
         Usage(argv[0]);
         return 0;
      } else if (!arg.empty() && arg[0] == '-') {
         std::cerr << "ramdump: unknown option '" << arg << "'\n";
         return 1;
      } else if (!file) {
         file = argv[i];
      } else {
         region = argv[i];
      }
   }

   if (!file) {
      Usage(argv[0]);
      return 1;
   }

   FILE *out = stdout;
   std::unique_ptr<FILE, int (*)(FILE *)> owned(nullptr, std::fclose);
   if (outPath) {
      out = std::fopen(outPath, "w");
      if (!out) {
         std::cerr << "ramdump: cannot write " << outPath << "\n";
         return 1;
      }
      owned.reset(out);
   }

   SamWriter writer(out);

   if (headerOnly) {
      WriteHeader(file, writer);
      return 0;
   }

   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   if (!reader) {
      std::cerr << "ramdump: cannot open " << file << "\n";
      return 1;
   }

   if (withHeader && !countOnly)
      WriteHeader(file, writer);

   // Counting with no filter never has to decode a record, so let the scan's own
   // count stand rather than materialising every row to throw it away.
   if (countOnly && requireFlags == 0 && excludeFlags == 0) {
      const Long64_t total = ramntuplescan(*reader, region, nullptr);
      writer.Int(total);
      writer.EndLine();
      return 0;
   }

   auto view = reader->GetView<RAMNTupleRecord>("record");
   Long64_t kept = 0;

   ramntuplescan(*reader, region, [&](Long64_t row) {
      const auto &rec = view(row);
      const uint16_t flag = rec.GetFLAG();
      if ((flag & requireFlags) != requireFlags)
         return;
      if (flag & excludeFlags)
         return;
      kept++;
      if (!countOnly)
         WriteRecord(rec, writer);
   });

   if (countOnly) {
      writer.Int(kept);
      writer.EndLine();
   }

   return 0;
}
