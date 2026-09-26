// Qualities written through QualityBlockWriter in small blocks must come back
// through QualityBlockReader exactly, whichever order the rows are read in.
#include "ramcore/QualityBlocks.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

struct Read {
   uint16_t flag;
   std::string qual;
};

const char *const kFile = "quality_blocks_test.root";

// Writes the reads in blocks of `block` records and returns the block sizes.
std::vector<uint32_t> Write(const std::vector<Read> &reads, std::size_t block, uint32_t policy)
{
   RAMNTupleRecord::InitializeRefs();
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   std::vector<uint32_t> sizes;
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      QualityBlockWriter out(
         writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
         [&writer](ROOT::REntry &e) { writer->Fill(e); }, block);
      for (const auto &read : reads) {
         RAMNTupleRecord &rec = out.Record();
         rec.SetBit(policy); // as the converters do: the bit is added, never replaced
         rec.SetQNAME("r");
         rec.SetFLAG(read.flag);
         rec.SetRNAME("chr1");
         rec.SetQUAL(read.qual);
         out.Add();
      }
      out.Finish();
      sizes = out.TakeBlockSizes();
   }
   std::vector<uint64_t> ends;
   uint64_t row = 0;
   for (const uint32_t n : sizes) {
      row += n;
      ends.push_back(row - 1);
   }
   RAMNTupleRecord::SetQualBlockEnds(ends);
   RAMNTupleRecord::WriteAllRefs(*file);
   file->Close();
   return sizes;
}

// Every read's quality, read back in the row order given.
void ExpectQualities(const std::vector<Read> &reads, const std::vector<std::size_t> &order)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   for (const std::size_t row : order)
      EXPECT_EQ(quals.Get(view(row), row), reads[row].qual) << "row " << row;
}

class QualityBlocksTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override { std::remove(kFile); }
};

TEST_F(QualityBlocksTest, EveryQualityComesBackInAnyOrder)
{
   const std::vector<Read> reads = {
      {0, "IIIIHHHGGF"},
      {0x10, "#####ABCDE"},      // reverse strand: modelled in sequencing order
      {0x80, "!~!~!~"},          // READ2, both ends of the SAM range
      {0, "*"},                  // no quality: stays in the record
      {0x90, "ABCDEFGHIJKLMNO"}, // reverse READ2, another length
      {0, "AB C"},               // a space is not SAM quality: stays in the record
      {0, "FFFFFFFFFF"},
      {0x10, "5"},
      {0, "*"},
      {0, "JJJJJJJJJJJJ"},
   };
   const auto sizes = Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   EXPECT_EQ(sizes, (std::vector<uint32_t>{3, 3, 3, 1}));

   std::vector<std::size_t> forward;
   std::vector<std::size_t> backward;
   for (std::size_t i = 0; i < reads.size(); i++) {
      forward.push_back(i);
      backward.push_back(reads.size() - 1 - i);
   }
   ExpectQualities(reads, forward);
   ExpectQualities(reads, backward);
}

TEST_F(QualityBlocksTest, OnlyPlainSamQualityGoesIntoTheBlock)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0, "*"}, {0, "AB C"}, {0, "IIII"}};
   Write(reads, /*block=*/10, RAMNTupleRecord::kPhred33);

   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   const std::vector<bool> inBlock = {true, false, false, true};
   for (std::size_t row = 0; row < reads.size(); row++) {
      const auto &rec = view(row);
      EXPECT_EQ(rec.TestBit(RAMNTupleRecord::kQualInBlock), inBlock[row]) << "row " << row;
      EXPECT_EQ(rec.qual.empty(), inBlock[row]) << "row " << row;
   }
}

TEST_F(QualityBlocksTest, ABlockWithoutQualitiesNeedsNoData)
{
   const std::vector<Read> reads = {{0, "*"}, {0, "*"}, {0, "*"}, {0, "HHHH"}};
   Write(reads, /*block=*/3, RAMNTupleRecord::kPhred33);
   ExpectQualities(reads, {3, 0, 1, 2});
}

TEST_F(QualityBlocksTest, LossyPoliciesKeepTheirOwnEncoding)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0, "!!!!"}};
   Write(reads, /*block=*/3, RAMNTupleRecord::kDrop);
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   QualityBlockReader quals(*reader);
   EXPECT_FALSE(view(0).TestBit(RAMNTupleRecord::kQualInBlock));
   EXPECT_EQ(quals.Get(view(0), 0), "*");
   EXPECT_EQ(quals.Get(view(1), 1), "*");
}

// The two record objects are reused; a record whose quality stays in the
// record must not keep the in-block bit of the record before it.
TEST_F(QualityBlocksTest, TheInBlockBitIsClearedOnReuse)
{
   const std::vector<Read> reads = {{0, "IIII"}, {0, "IIII"}, {0, "*"}, {0, "*"}};
   Write(reads, /*block=*/10, RAMNTupleRecord::kPhred33);
   ExpectQualities(reads, {0, 1, 2, 3});
}

} // namespace
