#ifndef RAMCORE_RAMNTUPLEVIEW_H
#define RAMCORE_RAMNTUPLEVIEW_H
#include <Rtypes.h>
struct RAMNTupleViewOpts {
   bool fCache = true;
   bool fPerfStats = false;
   std::string perfStatsFilename = "perf.root";
};
Long64_t ramntupleview(const char *file, const char *query = "", const RAMNTupleViewOpts & = RAMNTupleViewOpts());

#endif // RAMCORE_RAMNTUPLEVIEW_H
