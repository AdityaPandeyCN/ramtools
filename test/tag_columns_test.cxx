#include "ramcore/QualityBlocks.h"
#include "ramcore/TagColumns.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <TFile.h>

#include <gtest/gtest.h>

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace {

const char *const kFile = "tag_columns_test.root";

using Tags = std::vector<std::string>;

void Write(const std::vector<Tags> &records, const std::vector<TagColumn> &columns)
{
   RAMNTupleRecord::InitializeRefs();
   std::unique_ptr<TFile> file(TFile::Open(kFile, "RECREATE"));
   {
      auto model = RAMNTupleRecord::MakeModel();
      AddTagFields(*model, columns);
      auto writer = ROOT::RNTupleWriter::Append(std::move(model), "RAM", *file);
      QualityBlockWriter out(writer->GetModel().CreateEntry(), writer->GetModel().CreateEntry(),
                             [&writer](ROOT::REntry &e) { writer->Fill(e); });
      TagWriter tags(columns);
      for (const auto &t : records) {
         RAMNTupleRecord &rec = out.Record();
         rec.SetQNAME("r");
         rec.SetRNAME("chr1");
         rec.SetQUAL("*");
         rec.tags = t;
         tags.Move(rec, out.Entry());
         out.Add();
      }
      out.Finish();
      RAMNTupleRecord::SetQualBlockEnds(out.TakeBlockEnds());
   }
   RAMNTupleRecord::WriteAllRefs(*file);
   file->Close();
}

// Every record's tags and the ones left as text in the record.
void ExpectTags(const std::vector<Tags> &records, const std::vector<Tags> &text)
{
   auto reader = RAMNTupleRecord::OpenRAMFile(kFile);
   ASSERT_NE(reader, nullptr);
   auto view = reader->GetView<RAMNTupleRecord>("record");
   TagReader tags(*reader);
   for (std::size_t row = 0; row < records.size(); row++) {
      EXPECT_EQ(tags.Get(view(row), row), records[row]) << "row " << row;
      EXPECT_EQ(view(row).GetTags(), text[row]) << "row " << row;
   }
}

class TagColumnsTest : public ::testing::Test {
protected:
   // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
   void TearDown() override { std::remove(kFile); }
};

TEST_F(TagColumnsTest, SamplerKeepsTagsOfOnePercentOfRecordsInFirstSeenOrder)
{
   TagSampler sampler;
   for (int i = 0; i < 200; i++) {
      sampler.AddRecord();
      sampler.AddTag("NM:i:0");
      sampler.AddTag(i < 2 ? "XA:Z:a" : "MD:Z:10");
      sampler.AddTag(i < 1 ? "OC:Z:10M" : "XS:f:1.5"); // 0.5%; float has no column
   }
   const auto columns = sampler.Columns();
   ASSERT_EQ(columns.size(), 3U);
   EXPECT_EQ(columns[0].FieldName(), "tag_NM_i");
   EXPECT_EQ(columns[1].FieldName(), "tag_XA_Z");
   EXPECT_EQ(columns[2].FieldName(), "tag_MD_Z");
}

TEST_F(TagColumnsTest, TagsComeBackInTheirOrder)
{
   const std::vector<TagColumn> columns = {{"NM", 'i'}, {"MD", 'Z'}, {"XT", 'A'}, {"RG", 'Z'}};
   const std::vector<Tags> records = {
      {"NM:i:3", "MD:Z:0N0T48", "XT:A:R", "RG:Z:grp"},
      {"RG:Z:grp", "XS:f:1.5", "NM:i:-7", "BQ:Z:@@@"},
      {},
      {"NM:i:+5", "NM:i:05", "XT:A:RU", "OC:Z:10M"},
      {"NM:i:1", "NM:i:2"},
      {"NM:i:2147483647", "MD:Z:", "XB:B:c,1,2"},
      {"NM:i:2147483648"},
   };
   Write(records, columns);
   ExpectTags(records, {
                          {},
                          {"XS:f:1.5", "BQ:Z:@@@"},
                          {},
                          {"NM:i:+5", "NM:i:05", "XT:A:RU", "OC:Z:10M"}, // not as they print back
                          {"NM:i:2"},                                    // a key repeated in a record
                          {"XB:B:c,1,2"},
                          {"NM:i:2147483648"}, // beyond int32
                       });
}

TEST_F(TagColumnsTest, WithoutColumnsTagsStayInTheRecord)
{
   const std::vector<Tags> records = {{"NM:i:3", "RG:Z:grp"}, {}};
   Write(records, {});
   ExpectTags(records, records);
}

} // namespace
