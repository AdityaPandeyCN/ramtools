#ifndef RAMCORE_SAMTOTREE_H
#define RAMCORE_SAMTOTREE_H

#include <Rtypes.h> 

/// Converts a SAM file to the legacy TTree-backed RAM format.
///
/// \param compression_settings ROOT compression code, algorithm*100+level
///        (207 = LZMA-7, 505 = ZSTD-5, 404 = LZ4-4, 101 = ZLIB-1). This used to
///        be an algorithm id with the level pinned to 1, which made any
///        comparison against the RNTuple writer -- which takes the full code --
///        a comparison of two different levels.
void samtoram(const char *datafile,
              const char *treefile,
              bool index, bool split, bool cache,
              Int_t compression_settings,
              UInt_t quality_policy);

#endif // RAMCORE_SAMTOTREE_H 

