#ifndef RAMCORE_TAGCOLUMNS_H
#define RAMCORE_TAGCOLUMNS_H

// Common SAM tags in typed columns beside the record, one per key and type, so
// numbers are stored as numbers and each tag compresses on its own. The other
// tags stay as text in the record. "tagmask" says which columns a row has; a
// row whose tags are not in column order, text tags last, also has "taglayout".

#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RFieldToken.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

struct TagColumn {
   std::string key; ///< Two characters, such as "NM".
   char type = 'Z'; ///< 'i', 'A' or 'Z'.

   std::string FieldName() const { return "tag_" + key + "_" + type; }
};

/// Chooses the columns from the tags of the first records of an input.
class TagSampler {
public:
   void AddRecord() { m_records++; }
   void AddTag(std::string_view tag);
   /// Keys and types of type 'i', 'A' or 'Z' found in at least 1% of the records,
   /// in the order they first appear.
   std::vector<TagColumn> Columns() const;

private:
   std::map<std::string, std::size_t> m_counts; ///< By key and type, such as "NMi".
   std::vector<std::string> m_order;
   std::size_t m_records = 0;
};

/// Adds the mask, the layout and one field per column; nothing without columns.
void AddTagFields(ROOT::RNTupleModel &model, const std::vector<TagColumn> &columns);

/// Moves the tags that fit a column from the record to the entry's fields.
class TagWriter {
public:
   explicit TagWriter(std::vector<TagColumn> columns = {}) : m_columns(std::move(columns)) {}

   void Move(RAMNTupleRecord &rec, ROOT::REntry &entry);

private:
   std::vector<TagColumn> m_columns;
   std::vector<ROOT::RFieldToken> m_tokens; ///< Per column, then mask and layout; from the first entry.
   std::vector<bool> m_used;
};

/// Returns a record's tags as SAM text, in their original order.
class TagReader {
public:
   explicit TagReader(ROOT::RNTupleReader &reader);

   std::vector<std::string> Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   struct Column {
      TagColumn tag;
      std::optional<ROOT::RNTupleView<std::int32_t>> ints;
      std::optional<ROOT::RNTupleView<char>> chars;
      std::optional<ROOT::RNTupleView<std::string>> strings;
   };

   void Append(std::string &tag, Column &column, ROOT::NTupleSize_t row);

   std::optional<ROOT::RNTupleView<std::uint64_t>> m_mask;
   std::optional<ROOT::RNTupleView<std::string>> m_layout;
   std::vector<Column> m_columns;
};

#endif // RAMCORE_TAGCOLUMNS_H
