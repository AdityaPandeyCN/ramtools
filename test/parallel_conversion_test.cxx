// Many small blocks on several threads must equal one block on one thread.
#include <gtest/gtest.h>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>
#include <TFile.h>
#include <TList.h>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToNTuple.h"
#include "rntuple/RAMNTupleRecord.h"

namespace {

constexpr size_t kTinyBlock = 700;

std::string Header(const char *order = "coordinate")
{
   std::string h{};
   h += "@HD\tVN:1.6\tSO:";
   h += order;
   h += "\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n@SQ\tSN:chr3\tLN:100000\n";
   h += "@PG\tID:test\tPN:test\n";
   return h;
}

std::string Record(int n, const char *rname, int flag, int pos, const char *rnext = "*", int pnext = 0,
                   const char *cigar = "50M")
{
   const bool mapped = !(flag & 4);
   std::string r = "r" + std::to_string(n) + "\t" + std::to_string(flag) + "\t" + rname + "\t" + std::to_string(pos) +
                   "\t60\t" + (mapped ? cigar : "*") + "\t" + rnext + "\t" + std::to_string(pnext) + "\t0\t" +
                   (mapped ? std::string(50, 'A') : "*") + "\t" + (mapped ? std::string(50, 'F') : "*") +
                   "\tNM:i:" + std::to_string(n % 3) + "\tRG:Z:g" + std::to_string(n % 2) + "\n";
   return r;
}

std::string SortedRecords(int per_ref)
{
   std::string s{};
   int n = 0;
   for (int i = 0; i < per_ref; ++i)
      s += Record(n++, "chr1", (i % 2) ? 99 : 0, 1 + (i * 37), (i % 5) ? "=" : "chr2", 1 + (i * 41));
   s += Record(n++, "chr1", 4, 1 + (per_ref * 37) + 100); // placed unmapped read after the last mapped one
   for (int i = 0; i < per_ref; ++i)
      s += Record(n++, "chr2", 0, 5 + (i * 41), (i % 7) ? "=" : "chr1", 1 + (i * 37));
   s += Record(n++, "*", 4, 0);
   s += Record(n++, "*", 4, 0);
   return s;
}

void WriteSam(const char *path, const std::string &content)
{
   std::ofstream sam(path);
   sam << content;
}

struct Dump {
   std::vector<std::string> lines;
   bool sorted = false;
   uint32_t max_span = 0;
   std::vector<std::string> rname_refs;
};

Dump ReadBack(const char *file)
{
   Dump d;
   auto reader = RAMNTupleRecord::OpenRAMFile(file);
   EXPECT_NE(reader, nullptr) << file;
   if (!reader)
      return d;
   d.sorted = RAMNTupleRecord::IsCoordinateSorted();
   d.max_span = RAMNTupleRecord::GetMaxRefSpan();
   d.rname_refs = RAMNTupleRecord::GetRnameRefs()->GetRefs();
   auto view = reader->GetView<RAMNTupleRecord>("record");
   for (auto i : reader->GetEntryRange()) {
      const RAMNTupleRecord &r = view(i);
      std::string line = r.GetQNAME() + "\t" + std::to_string(r.GetFLAG()) + "\t" + r.GetRNAME() + "\t" +
                         std::to_string(r.GetPOS()) + "\t" + std::to_string(r.GetMAPQ()) + "\t" + r.GetCIGAR() +
                         "\t" + r.GetRNEXT() + "\t" + std::to_string(r.GetPNEXT()) + "\t" +
                         std::to_string(r.GetTLEN()) + "\t" + r.GetSEQ() + "\t" + r.GetQUAL();
      for (const auto &t : r.GetTags())
         line += "\t" + t;
      d.lines.push_back(line);
   }
   return d;
}

class ParallelConversionTest : public ::testing::Test {
protected:
   static constexpr const char *kSam = "parallel_test.sam";
   static constexpr const char *kSeq = "parallel_test_seq.root";
   static constexpr const char *kPar = "parallel_test_par.root";

   void TearDown() override
   {
      std::remove(kSam);
      std::remove(kSeq);
      std::remove(kPar);
   }
};

TEST_F(ParallelConversionTest, ManyBlocksOnSeveralThreadsMatchOneBlockRecordForRecord)
{
   WriteSam(kSam, Header() + SortedRecords(300));
   samtoramntuple(kSam, kSeq, 505, 0);
   const Dump seq = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, 505, 0, 3, kTinyBlock);
   const Dump par = ReadBack(kPar);

   ASSERT_EQ(par.lines.size(), 603U);
   EXPECT_EQ(par.lines, seq.lines) << "same records in the same order as with one block";
   EXPECT_TRUE(par.sorted);
   EXPECT_EQ(par.sorted, seq.sorted);
   EXPECT_EQ(par.max_span, seq.max_span);
   EXPECT_EQ(par.max_span, 50U);
   EXPECT_EQ(par.rname_refs, seq.rname_refs) << "reference ids follow the header in both";

   auto reader = ROOT::RNTupleReader::Open("RAM", kPar);
   EXPECT_GT(reader->GetDescriptor().GetNClusters(), 20U);
   const RAMNTupleViewOpts opts = {true, false, ""};
   EXPECT_EQ(ramntupleview(kPar, "chr1:11150-11250", opts), ramntupleview(kSeq, "chr1:11150-11250", opts));
   EXPECT_EQ(ramntupleview(kPar, "chr1:11150-11250", opts), 1) << "the placed unmapped read";
   EXPECT_EQ(ramntupleview(kPar, "chr2:4000-8200", opts), ramntupleview(kSeq, "chr2:4000-8200", opts));
   EXPECT_EQ(ramntupleview(kPar, "chr2", opts), 300);
   EXPECT_EQ(ramntupleview(kPar, "chr3", opts), 0);
}

TEST_F(ParallelConversionTest, OneThreadAndOneBlockGiveOneCluster)
{
   WriteSam(kSam, Header() + SortedRecords(40));
   samtoramntuple(kSam, kSeq, 505, 0);
   const Dump seq = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, 505, 0, 1);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_TRUE(par.sorted);

   auto reader = ROOT::RNTupleReader::Open("RAM", kPar);
   EXPECT_EQ(reader->GetDescriptor().GetNClusters(), 1U);
}

// Each block is sorted on its own; only the step between blocks goes back.
TEST_F(ParallelConversionTest, DisorderAcrossABlockBoundaryMarksTheFileUnsorted)
{
   std::string s = Header("unsorted");
   int n = 0;
   for (int i = 0; i < 12; ++i)
      s += Record(n++, "chr1", 0, 5000 + (i * 10)); // 12 records at ~150 bytes each: a few blocks
   for (int i = 0; i < 12; ++i)
      s += Record(n++, "chr1", 0, 1000 + (i * 10)); // sorted within itself, but before the block above
   WriteSam(kSam, s);

   testing::internal::CaptureStderr();
   samtoramntuple(kSam, kPar, 505, 0, 2, kTinyBlock);
   EXPECT_NE(testing::internal::GetCapturedStderr().find("not in coordinate order"), std::string::npos);

   const Dump par = ReadBack(kPar);
   EXPECT_FALSE(par.sorted);
   ASSERT_EQ(par.lines.size(), 24U);
   EXPECT_EQ(par.lines.front().substr(0, 3), "r0\t");
   EXPECT_EQ(par.lines.back().substr(0, 4), "r23\t");
   const RAMNTupleViewOpts opts = {true, false, ""};
   EXPECT_EQ(ramntupleview(kPar, "chr1:1000-1200", opts), 12);
}

TEST_F(ParallelConversionTest, UnplacedRecordsBeforePlacedOnesInALaterBlockMeanUnsorted)
{
   std::string s = Header("unsorted");
   int n = 0;
   for (int i = 0; i < 10; ++i)
      s += Record(n++, "*", 4, 0);
   for (int i = 0; i < 10; ++i)
      s += Record(n++, "chr1", 0, 1000 + (i * 10));
   WriteSam(kSam, s);

   testing::internal::CaptureStderr();
   samtoramntuple(kSam, kPar, 505, 0, 2, kTinyBlock);
   testing::internal::GetCapturedStderr();

   const Dump par = ReadBack(kPar);
   EXPECT_FALSE(par.sorted);
   EXPECT_EQ(par.lines.size(), 20U);
}

TEST_F(ParallelConversionTest, HeaderLongerThanABlockIsReadBeforeTheRecords)
{
   std::string s = Header();
   for (int i = 0; i < 40; ++i)
      s += "@CO\tcomment line number " + std::to_string(i) + " padding padding padding padding\n";
   s += SortedRecords(20);
   WriteSam(kSam, s);
   samtoramntuple(kSam, kSeq, 505, 0);
   const Dump seq = ReadBack(kSeq);

   samtoramntuple(kSam, kPar, 505, 0, 2, kTinyBlock);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_EQ(par.rname_refs, seq.rname_refs);

   auto file = std::unique_ptr<TFile>(TFile::Open(kPar));
   auto *headers = file->Get<TList>("headers");
   ASSERT_NE(headers, nullptr);
   EXPECT_EQ(headers->GetSize(), 45) << "@HD, three @SQ, @PG and forty @CO lines";
}

TEST_F(ParallelConversionTest, MalformedAndEmptyLinesAreSkippedInEveryBlock)
{
   std::string s = Header();
   s += Record(0, "chr1", 0, 100);
   s += "\n";
   s += "bad\tx\tchr1\t200\t60\t50M\t*\t0\t0\t*\t*\n"; // flag is not a number
   s += Record(1, "chr1", 0, 300);
   s += "short\t0\tchr1\n";
   s += Record(2, "chr1", 0, 400, "*", 0, "10M5Q35M"); // bad CIGAR operator
   s += Record(3, "chr1", 0, 500);
   s += "@CO\ta comment among the records\n";
   s += Record(4, "chr1", 0, 600);
   WriteSam(kSam, s);

   testing::internal::CaptureStderr();
   samtoramntuple(kSam, kSeq, 505, 0);
   samtoramntuple(kSam, kPar, 505, 0, 2, kTinyBlock);
   testing::internal::GetCapturedStderr();

   const Dump seq = ReadBack(kSeq);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(seq.lines.size(), 4U);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_TRUE(par.sorted);

   auto file = std::unique_ptr<TFile>(TFile::Open(kPar));
   auto *headers = file->Get<TList>("headers");
   ASSERT_NE(headers, nullptr);
   EXPECT_EQ(headers->GetSize(), 6) << "the late @CO line is kept";
}

TEST_F(ParallelConversionTest, LastLineWithoutNewlineIsKept)
{
   std::string s = Header() + SortedRecords(30);
   s.pop_back();
   WriteSam(kSam, s);
   samtoramntuple(kSam, kSeq, 505, 0);
   const Dump seq = ReadBack(kSeq);
   ASSERT_EQ(seq.lines.size(), 63U);

   samtoramntuple(kSam, kPar, 505, 0, 2, kTinyBlock);
   const Dump par = ReadBack(kPar);
   EXPECT_EQ(par.lines, seq.lines);
   EXPECT_TRUE(par.sorted);
}

// A directory opens but fails on read.
TEST_F(ParallelConversionTest, ReadErrorIsReportedAfterTheWorkersAreStopped)
{
   testing::internal::CaptureStdout();
   EXPECT_FALSE(samtoramntuple(".", kPar, 505, 0, 3, kTinyBlock));
   EXPECT_NE(testing::internal::GetCapturedStdout().find("Failed to read SAM file"), std::string::npos);
}

TEST_F(ParallelConversionTest, MissingInputIsReported)
{
   testing::internal::CaptureStdout();
   EXPECT_FALSE(samtoramntuple("does_not_exist.sam", kPar, 505, 0, 2));
   EXPECT_NE(testing::internal::GetCapturedStdout().find("Failed to parse SAM file"), std::string::npos);
}

} // namespace
