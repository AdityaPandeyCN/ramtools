// raminfo: describe a RAM file and show where its bytes go, one row per
// record member, and with -v one row per physical column.

#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RColumnElementBase.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <ROOT/RNTupleInspector.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <TFile.h>
#include <TList.h>

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

using ROOT::Experimental::RNTupleInspector;

double MB(std::uint64_t bytes)
{
   return static_cast<double>(bytes) / 1e6;
}

double Ratio(std::uint64_t uncompressed, std::uint64_t compressed)
{
   return compressed > 0 ? static_cast<double>(uncompressed) / static_cast<double>(compressed) : 0.0;
}

void PrintSummary(const char *file, const RNTupleInspector &inspector)
{
   const auto &desc = inspector.GetDescriptor();
   std::unique_ptr<TFile> f(TFile::Open(file, "READ"));
   const auto fileBytes = f ? static_cast<std::uint64_t>(f->GetSize()) : 0;
   int headerLines = -1;
   if (f) {
      if (auto *headers = f->Get<TList>("headers"))
         headerLines = headers->GetSize();
   }

   printf("file          %s, %.2f MB\n", file, MB(fileBytes));
   printf("records       %llu in %zu clusters\n", static_cast<unsigned long long>(desc.GetNEntries()),
          static_cast<std::size_t>(desc.GetNClusters()));
   printf("compression   %s\n", inspector.GetCompressionSettingsAsString().c_str());
   printf("RAM ntuple    %.2f MB compressed, %.2f MB uncompressed, ratio %.2f\n", MB(inspector.GetCompressedSize()),
          MB(inspector.GetUncompressedSize()), Ratio(inspector.GetUncompressedSize(), inspector.GetCompressedSize()));
   printf("coordinate    %s; longest reference span %u\n",
          RAMNTupleRecord::IsCoordinateSorted() ? "sorted" : "not sorted", RAMNTupleRecord::GetMaxRefSpan());
   printf("references    %zu RNAME, %zu RNEXT\n", RAMNTupleRecord::GetRnameRefs()->Size(),
          RAMNTupleRecord::GetRnextRefs()->Size());
   if (headerLines >= 0)
      printf("header        %d lines\n", headerLines);
   printf("\n");
}

void PrintColumns(const RNTupleInspector &inspector, ROOT::DescriptorId_t fieldId, std::uint64_t total)
{
   const auto &desc = inspector.GetDescriptor();
   for (const auto &column : desc.GetColumnIterable(fieldId)) {
      if (column.IsAliasColumn())
         continue;
      const auto &ci = inspector.GetColumnInspector(column.GetPhysicalId());
      const std::string name = "  " + desc.GetQualifiedFieldName(fieldId) + " #" + std::to_string(column.GetIndex());
      printf("%-30s %12.2f %6.1f%% %12.2f %7.2f %9llu  %s\n", name.c_str(), MB(ci.GetCompressedSize()),
             100.0 * static_cast<double>(ci.GetCompressedSize()) / static_cast<double>(total),
             MB(ci.GetUncompressedSize()), Ratio(ci.GetUncompressedSize(), ci.GetCompressedSize()),
             static_cast<unsigned long long>(ci.GetNPages()),
             ROOT::Internal::RColumnElementBase::GetColumnTypeName(ci.GetType()));
   }
   for (const auto &child : desc.GetFieldIterable(fieldId))
      PrintColumns(inspector, child.GetId(), total);
}

void PrintMembers(const RNTupleInspector &inspector, bool verbose)
{
   const auto &desc = inspector.GetDescriptor();
   const auto recordId = desc.FindFieldId("record");
   const auto total = inspector.GetCompressedSize();
   const auto entries = desc.GetNEntries();

   printf("%-30s %12s %7s %12s %7s %9s\n", "member", "MB", "share", "raw MB", "ratio", "B/record");
   if (verbose)
      printf("%-30s %12s %7s %12s %7s %9s  %s\n", "  column", "", "", "", "", "pages", "type");
   for (const auto &member : desc.GetFieldIterable(recordId)) {
      const auto &fi = inspector.GetFieldTreeInspector(member.GetId());
      printf("%-30s %12.2f %6.1f%% %12.2f %7.2f %9.2f\n", desc.GetQualifiedFieldName(member.GetId()).c_str(),
             MB(fi.GetCompressedSize()),
             100.0 * static_cast<double>(fi.GetCompressedSize()) / static_cast<double>(total),
             MB(fi.GetUncompressedSize()), Ratio(fi.GetUncompressedSize(), fi.GetCompressedSize()),
             entries > 0 ? static_cast<double>(fi.GetCompressedSize()) / static_cast<double>(entries) : 0.0);
      if (verbose)
         PrintColumns(inspector, member.GetId(), total);
   }
}

} // namespace

int main(int argc, char *argv[])
{
   const std::vector<std::string> args(argv + 1, argv + argc);
   bool verbose = false;
   std::string file{};
   bool usage = false;
   for (const auto &arg : args) {
      if (arg == "-v")
         verbose = true;
      else if (arg.empty() || arg[0] == '-' || !file.empty())
         usage = true;
      else
         file = arg;
   }
   if (usage || file.empty()) {
      fprintf(stderr, "Usage: raminfo [-v] <file.ram>\n"
                      "Prints the file's layout and the compressed size of each record member.\n"
                      "  -v   also list every column: its size, pages and on-disk type\n");
      return 1;
   }

   // Loads the name tables, the sort flag and the longest span.
   if (!RAMNTupleRecord::OpenRAMFile(file))
      return 1;
   auto inspector = RNTupleInspector::Create("RAM", file);

   PrintSummary(file.c_str(), *inspector);
   PrintMembers(*inspector, verbose);
   return 0;
}
