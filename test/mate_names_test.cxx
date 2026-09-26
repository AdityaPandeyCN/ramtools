#include "ramcore/MateNames.h"
#include "ramcore/QualityBlocks.h"
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
#include <string>
#include <vector>

namespace {

const char *const kFile = "mate_names_test.root";

// Writes the names, calling Reset() before the rows listed in `resets`.
void Write(const std::vector<std::string> &names, const std::vector<std::size_t> &resets = {})
{
   RAMNTupleRecord::InitializeRefs();
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   {
      auto writer = ROOT::RNTupleWriter::Append(RAMNTupleRecord::MakeModel(), "RAM", *file);
      QualityBlockWriter out(writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
                             [&writer](ROOT::REntry &e) { writer->Fill(e); });
      MateNameWriter mates;
      for (std::size_t row = 0; row < names.size(); row++) {
         if (std::find(resets.begin(), resets.end(), row) != resets.end())
            mates.Reset();
         RAMNTupleRecord &rec = out.Record();
         rec.SetQNAME(names[row]);
         rec.SetRNAME("chr1");
         rec.SetQUAL("*");
         mates.Move(rec, out.Entry());
         out.Add();
      }
      out.Finish();
      RAMNTupleRecord::SetQualBlockEnds(out.TakeBlockEnds());
   }
   RAMNTupleRecord::WriteAllRefs(*file);
   file->Close();
}

// Checks every name, and returns the names stored as text in the records.
std::vector<std::string> ReadBack(const std::vector<std::string> &names)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   EXPECT_NE(reader, nullptr);
   if (!reader)
      return {};
   auto view = reader->GetView<RAMNTupleRecord>("record");
   MateNameReader mates(*reader);
   std::vector<std::string> stored;
   for (std::size_t row = 0; row < names.size(); row++) {
      EXPECT_EQ(mates.Get(view(row), row), names[row]) << "row " << row;
      stored.push_back(view(row).GetQNAME());
   }
   return stored;
}

class MateNamesTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override { std::remove(kFile); }
};

TEST_F(MateNamesTest, MatesReferToTheFirstRead)
{
   const std::vector<std::string> names = {"a", "b", "a", "c", "b", "a", "a"};
   Write(names);
   // The third "a" starts a new pair with the fourth.
   EXPECT_EQ(ReadBack(names), (std::vector<std::string>{"a", "b", "", "c", "", "a", ""}));
}

TEST_F(MateNamesTest, NamesAreLookedForInTheLastWindowOnly)
{
   std::vector<std::string> names;
   for (std::size_t i = 0; i < 2 * kMateNameWindow + 10; i++)
      names.push_back("r" + std::to_string(i));
   names[kMateNameWindow + 5] = "r1";     // in the window before the current one
   names[2 * kMateNameWindow + 5] = "r2"; // two windows on: too far back
   Write(names);
   const auto stored = ReadBack(names);
   EXPECT_EQ(stored[kMateNameWindow + 5], "");
   EXPECT_EQ(stored[2 * kMateNameWindow + 5], "r2");
}

TEST_F(MateNamesTest, ResetForgetsEarlierRows)
{
   const std::vector<std::string> names = {"a", "b", "a", "b"};
   Write(names, /*resets=*/{2});
   EXPECT_EQ(ReadBack(names), names);
}

} // namespace
