#include "ramcore/SamParser.h"
#include <cctype>
#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>

namespace ramcore {

namespace {

bool ParseInt(const char *value, int &out, const char *field_name, size_t line_number, int min_value, int max_value)
{
   if (!value || *value == '\0') {
      std::cerr << "[SamParser] Warning: empty value for field '" << field_name << "' at line " << line_number
                << "; record skipped.\n";
      return false;
   }
   if (std::isspace(static_cast<unsigned char>(*value))) {
      std::cerr << "[SamParser] Warning: invalid integer value '" << value << "' for field '" << field_name
                << "' at line " << line_number << "; record skipped.\n";
      return false;
   }

   errno = 0;
   char *end = nullptr;
   long parsed = std::strtol(value, &end, 10);

   if (end == value || *end != '\0') {
      std::cerr << "[SamParser] Warning: invalid integer value '" << value << "' for field '" << field_name
                << "' at line " << line_number << "; record skipped.\n";
      return false;
   }

   if (errno == ERANGE || parsed < min_value || parsed > max_value) {
      std::cerr << "[SamParser] Warning: integer value '" << value << "' for field '" << field_name
                << "' is out of range at line " << line_number << "; record skipped.\n";
      return false;
   }

   out = static_cast<int>(parsed);
   return true;
}

// Every mandatory column holds at least one character, so an empty one means the
// line is malformed -- usually a stray tab shifting the later columns.
bool RequireNonEmpty(const char *value, const char *field_name, size_t line_number)
{
   if (value && *value != '\0') {
      return true;
   }
   std::cerr << "[SamParser] Warning: empty value for field '" << field_name << "' at line " << line_number
             << "; record skipped.\n";
   return false;
}

} // namespace

// NOLINTNEXTLINE(misc-use-internal-linkage)
void StripCRLF(char *str)
{
   size_t len = strlen(str);
   // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
   while (len > 0 && (str[len - 1] == '\n' || str[len - 1] == '\r')) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      str[--len] = '\0';
   }
}

bool SamParser::ParseFile(const char *filename, HeaderCallback header_cb, RecordCallback record_cb)
{
   FILE *fp = fopen(filename, "r");
   if (!fp) {
      return false;
   }

   lines_processed_ = 0;
   records_processed_ = 0;

   // getline() grows its buffer to fit the line. A fixed buffer splits longer
   // records into fragments that then parse as separate malformed records, so
   // the read is lost outright -- and PacBio/ONT lines routinely exceed 100 kB.
   char *line = nullptr;
   size_t capacity = 0;

   SamRecord record;

   while (getline(&line, &capacity, fp) != -1) {
      lines_processed_++;

      StripCRLF(line);

      if (line[0] == '\0') {
         continue;
      }

      if (line[0] == '@') {
         char *tab = strchr(line, '\t');
         if (tab) {
            *tab = '\0';
            if (header_cb) {
               header_cb(line, tab + 1);
            }
         } else {
            if (header_cb) {
               header_cb(line, "");
            }
         }
         continue;
      }

      record.Clear();
      if (ParseLine(line, record)) {
         if (record_cb) {
            record_cb(record, records_processed_);
         }
         records_processed_++;
      }
   }

   free(line);
   fclose(fp);
   return true;
}

bool SamParser::ParseLine(char *line, SamRecord &record)
{
   int field_num = 0;
   char *cursor = line;
   bool last_field = false;
   static const char *kMandatory[] = {"qname", "flag", "rname", "pos",  "mapq", "cigar",
                                      "rnext", "pnext", "tlen", "seq", "qual"};

   // Split on TAB without collapsing runs of delimiters: "a\t\tb" is three
   // fields, the middle one empty. strtok() reported two and shifted every later
   // column left, so QUAL would be read as SEQ.
   while (!last_field) {
      char *token = cursor;
      char *tab = strchr(cursor, '\t');
      if (tab) {
         *tab = '\0';
         cursor = tab + 1;
      } else {
         last_field = true;
      }

      if (field_num <= 10 && !RequireNonEmpty(token, kMandatory[field_num], lines_processed_)) {
         return false;
      }

      switch (field_num) {
      case 0: record.qname = token; break;
      case 1:
         if (!ParseInt(token, record.flag, "flag", lines_processed_, 0, std::numeric_limits<unsigned short>::max()))
            return false;
         break;
      case 2: record.rname = token; break;
      case 3:
         if (!ParseInt(token, record.pos, "pos", lines_processed_, 0, std::numeric_limits<int>::max()))
            return false;
         break;
      case 4:
         if (!ParseInt(token, record.mapq, "mapq", lines_processed_, 0, std::numeric_limits<unsigned char>::max()))
            return false;
         break;
      case 5: record.cigar = token; break;
      case 6: record.rnext = token; break;
      case 7:
         if (!ParseInt(token, record.pnext, "pnext", lines_processed_, 0, std::numeric_limits<int>::max()))
            return false;
         break;
      case 8:
         if (!ParseInt(token, record.tlen, "tlen", lines_processed_, std::numeric_limits<int>::min() + 1,
                       std::numeric_limits<int>::max()))
            return false;
         break;
      case 9: record.seq = token; break;
      case 10: record.qual = token; break;
      default: record.optional_fields.push_back(token); break;
      }

      field_num++;
   }

   if (field_num < 11) {
      std::cerr << "[SamParser] Warning: line " << lines_processed_ << " has " << field_num
                << " fields, at least 11 are required; record skipped.\n";
      return false;
   }

   return true;
}

} // namespace ramcore
