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

/// Calls on_row for every record overlapping query, in file order, and returns
/// how many there were. The overlap rule is samtools': same reference, and
/// [POS, POS+refspan-1] meets the query interval. An empty query or "*" visits
/// the whole file. The reader must come from RAMNTupleRecord::OpenRAMFile so
/// the reference names and the index are loaded; on_row may be null.
Long64_t ramntuplescan(ROOT::RNTupleReader &reader, const char *query, const std::function<void(Long64_t)> &on_row);

/// Counts the records overlapping query.
Long64_t ramntupleview(const char *file, const char *query = "", const RAMNTupleViewOpts & = RAMNTupleViewOpts());

#endif // RAMCORE_RAMNTUPLEVIEW_H
