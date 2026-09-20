#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstdint>
#include <string>

/// Converts a SAM file to one RAM file on \p threads worker threads. The input
/// is read in blocks of \p block_bytes whole lines; each worker parses, encodes
/// and compresses its own blocks through ROOT's RNTupleParallelWriter. Records
/// keep their input order: every block becomes one or more clusters and the
/// clusters are appended in block order, so a coordinate-sorted input gives a
/// file that region queries can seek in. The default block size is a good
/// value; tests use a small one to exercise many blocks.
void samtoramntuple(const char *datafile, const char *treefile, int compression_algorithm, uint32_t quality_policy,
                    int threads = 1, size_t block_bytes = 64u << 20);

/// Writes one RAM file per reference, named <output_prefix>_<rname>.root.
/// Records stream to their file as they are parsed, in the order they arrive.
void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy);

#endif
