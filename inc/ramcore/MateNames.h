#ifndef RAMCORE_MATENAMES_H
#define RAMCORE_MATENAMES_H

// A read's name is usually repeated a few rows later by its mate. The later
// record stores an empty name and, in "qnameref", how many rows back the row
// with the name is.

#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RFieldToken.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>

/// Rows a name is looked for in: at least this many, at most twice as many.
constexpr std::size_t kMateNameWindow = 100000;

/// Replaces a record's name by a reference when a recent row has the same name.
class MateNameWriter {
public:
   void Move(RAMNTupleRecord &rec, ROOT::REntry &entry);
   /// Forgets the rows so far; the next row does not follow them in the file.
   void Reset();

private:
   std::optional<ROOT::RFieldToken> m_token;
   std::unordered_map<std::string, std::uint64_t> m_recent; ///< Name -> row, since m_recentStart.
   std::unordered_map<std::string, std::uint64_t> m_older;  ///< The window before.
   std::uint64_t m_row = 0;
   std::uint64_t m_recentStart = 0;
};

/// Returns QNAME for records read from a RAM file.
class MateNameReader {
public:
   explicit MateNameReader(ROOT::RNTupleReader &reader);

   std::string Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   std::optional<ROOT::RNTupleView<std::uint32_t>> m_ref;
   std::optional<ROOT::RNTupleView<std::string>> m_names;
};

#endif // RAMCORE_MATENAMES_H
