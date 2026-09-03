#ifndef RAMCORE_RAMNTUPLEVIEW_H
#define RAMCORE_RAMNTUPLEVIEW_H

#include <Rtypes.h>

#include <functional>
#include <string>

namespace ROOT {
class RNTupleReader;
}

struct RAMNTupleViewOpts {
   bool fCache = true;
   bool fPerfStats = false;
   std::string perfStatsFilename = "perf.root";
};

/// Visits every row overlapping \a query in reference order and returns how many
/// there were.
///
/// The overlap rule is samtools': a record is in the region when its reference
/// matches and [POS, POS+refspan-1] intersects the query interval, secondary,
/// supplementary and placed-unmapped records included. An empty query or "*"
/// visits the whole file.
///
/// \a reader must come from RAMNTupleRecord::OpenRAMFile so that the reference
/// names and the index are loaded. \a on_row receives each matching row number;
/// the caller owns whatever view it reads them through.
Long64_t ramntuplescan(ROOT::RNTupleReader &reader, const char *query, const std::function<void(Long64_t)> &on_row);

/// Counts the records overlapping \a query -- ramntuplescan with a counting callback.
Long64_t ramntupleview(const char *file, const char *query = "", const RAMNTupleViewOpts & = RAMNTupleViewOpts());

#endif // RAMCORE_RAMNTUPLEVIEW_H
