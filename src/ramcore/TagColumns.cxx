#include "ramcore/TagColumns.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>

#include <charconv>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace {

constexpr const char *kMaskField = "tagmask";
constexpr const char *kLayoutField = "taglayout";
constexpr std::size_t kMaxColumns = 64;
// A layout entry for a tag kept as text in the record; tag keys start with a letter.
constexpr char kTextTag = '*';

bool HasColumnType(char type)
{
   return type == 'i' || type == 'A' || type == 'Z';
}

// "NM:i:0" -> key "NM", type 'i', value "0"; false when the tag is not of that form.
bool SplitTag(std::string_view tag, std::string_view &key, char &type, std::string_view &value)
{
   if (tag.size() < 5 || tag[2] != ':' || tag[4] != ':')
      return false;
   key = tag.substr(0, 2);
   type = tag[3];
   value = tag.substr(5);
   return true;
}

// Only integers that print back the same way go into a column, so "+5" or "05" stay text.
bool ParseInt(std::string_view value, std::int32_t &out)
{
   const char *end = value.data() + value.size();
   const auto [ptr, ec] = std::from_chars(value.data(), end, out);
   return ec == std::errc() && ptr == end && std::to_string(out) == value;
}

} // namespace

void TagSampler::AddTag(std::string_view tag)
{
   std::string_view key;
   std::string_view value;
   char type = 0;
   if (!SplitTag(tag, key, type, value) || !HasColumnType(type))
      return;
   std::string keyType(key);
   keyType.push_back(type);
   if (m_counts[keyType]++ == 0)
      m_order.push_back(keyType);
}

std::vector<TagColumn> TagSampler::Columns() const
{
   std::vector<TagColumn> columns;
   for (const auto &keyType : m_order) {
      if (m_counts.at(keyType) * 100 >= m_records && columns.size() < kMaxColumns)
         columns.push_back({keyType.substr(0, 2), keyType[2]});
   }
   return columns;
}

void AddTagFields(ROOT::RNTupleModel &model, const std::vector<TagColumn> &columns)
{
   if (columns.empty())
      return;
   model.MakeField<std::uint64_t>(kMaskField);
   model.MakeField<std::string>(kLayoutField);
   for (const auto &c : columns) {
      switch (c.type) {
      case 'i': model.MakeField<std::int32_t>(c.FieldName()); break;
      case 'A': model.MakeField<char>(c.FieldName()); break;
      default: model.MakeField<std::string>(c.FieldName()); break;
      }
   }
}

void TagWriter::Move(RAMNTupleRecord &rec, ROOT::REntry &entry)
{
   if (m_columns.empty())
      return;
   if (m_tokens.empty()) {
      for (const auto &c : m_columns)
         m_tokens.push_back(entry.GetToken(c.FieldName()));
      m_tokens.push_back(entry.GetToken(kMaskField));
      m_tokens.push_back(entry.GetToken(kLayoutField));
   }

   for (std::size_t c = 0; c < m_columns.size(); c++) {
      switch (m_columns[c].type) {
      case 'i': *entry.GetPtr<std::int32_t>(m_tokens[c]) = 0; break;
      case 'A': *entry.GetPtr<char>(m_tokens[c]) = 0; break;
      default: entry.GetPtr<std::string>(m_tokens[c])->clear(); break;
      }
   }
   std::uint64_t &mask = *entry.GetPtr<std::uint64_t>(m_tokens[m_columns.size()]);
   std::string &layout = *entry.GetPtr<std::string>(m_tokens.back());
   mask = 0;
   layout.clear();
   m_used.assign(m_columns.size(), false);

   // In column order with any text tags last, the mask alone gives the order.
   bool inOrder = true;
   bool hadText = false;
   std::size_t next = 0;
   std::size_t kept = 0;
   for (std::size_t t = 0; t < rec.tags.size(); t++) {
      std::string_view key;
      std::string_view value;
      char type = 0;
      bool moved = false;
      if (SplitTag(rec.tags[t], key, type, value)) {
         for (std::size_t c = 0; c < m_columns.size() && !moved; c++) {
            // A key repeated in one record keeps its later copies as text.
            if (m_used[c] || m_columns[c].type != type || m_columns[c].key != key)
               continue;
            switch (type) {
            case 'i': {
               std::int32_t &v = *entry.GetPtr<std::int32_t>(m_tokens[c]);
               moved = ParseInt(value, v);
               if (!moved)
                  v = 0;
               break;
            }
            case 'A':
               moved = value.size() == 1;
               if (moved)
                  *entry.GetPtr<char>(m_tokens[c]) = value[0];
               break;
            default:
               entry.GetPtr<std::string>(m_tokens[c])->assign(value);
               moved = true;
               break;
            }
            if (moved) {
               m_used[c] = true;
               mask |= std::uint64_t{1} << c;
               layout.append(key).push_back(type);
               inOrder = inOrder && !hadText && c >= next;
               next = c + 1;
            }
         }
      }
      if (!moved) {
         layout.push_back(kTextTag);
         hadText = true;
         if (kept != t)
            rec.tags[kept] = std::move(rec.tags[t]);
         kept++;
      }
   }
   rec.tags.resize(kept);
   if (inOrder)
      layout.clear();
}

TagReader::TagReader(ROOT::RNTupleReader &reader)
{
   const auto &desc = reader.GetDescriptor();
   if (desc.FindFieldId(kLayoutField) == ROOT::kInvalidDescriptorId)
      return;
   m_mask.emplace(reader.GetView<std::uint64_t>(kMaskField));
   m_layout.emplace(reader.GetView<std::string>(kLayoutField));
   for (const auto &field : desc.GetTopLevelFields()) {
      const std::string &name = field.GetFieldName();
      if (name.size() != 8 || name.rfind("tag_", 0) != 0 || name[6] != '_')
         continue;
      Column c;
      c.tag = {name.substr(4, 2), name[7]};
      switch (c.tag.type) {
      case 'i': c.ints.emplace(reader.GetView<std::int32_t>(name)); break;
      case 'A': c.chars.emplace(reader.GetView<char>(name)); break;
      default: c.strings.emplace(reader.GetView<std::string>(name)); break;
      }
      m_columns.push_back(std::move(c));
   }
}

void TagReader::Append(std::string &tag, Column &column, ROOT::NTupleSize_t row)
{
   tag.append(column.tag.key).push_back(':');
   tag.push_back(column.tag.type);
   tag.push_back(':');
   if (column.ints)
      tag += std::to_string((*column.ints)(row));
   else if (column.chars)
      tag.push_back((*column.chars)(row));
   else
      tag += (*column.strings)(row);
}

std::vector<std::string> TagReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!m_layout)
      return rec.GetTags();

   std::vector<std::string> tags;
   const std::string &layout = (*m_layout)(row);
   if (layout.empty()) {
      const std::uint64_t mask = (*m_mask)(row);
      for (std::size_t c = 0; c < m_columns.size(); c++) {
         if (mask & (std::uint64_t{1} << c))
            Append(tags.emplace_back(), m_columns[c], row);
      }
      tags.insert(tags.end(), rec.GetTags().begin(), rec.GetTags().end());
      return tags;
   }

   std::size_t text = 0;
   for (std::size_t i = 0; i < layout.size();) {
      if (layout[i] == kTextTag) {
         tags.push_back(rec.GetTags().at(text++));
         i++;
         continue;
      }
      const std::string_view key(layout.data() + i, 2);
      const char type = layout.at(i + 2);
      i += 3;
      Column *column = nullptr;
      for (auto &c : m_columns) {
         if (c.tag.type == type && c.tag.key == key)
            column = &c;
      }
      if (!column)
         throw std::runtime_error("record " + std::to_string(row) + " has a tag with no column: " + std::string(key));
      Append(tags.emplace_back(), *column, row);
   }
   return tags;
}
