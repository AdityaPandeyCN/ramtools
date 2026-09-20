#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstdint>
#include <string>

void samtoramntuple(const char *datafile, const char *treefile, bool split, bool cache, int compression_algorithm,
                    uint32_t quality_policy);

/// Converts with \p threads worker threads, each parsing, encoding and
/// compressing its own blocks of the input through ROOT's RNTupleParallelWriter.
/// Records keep their input order: every block becomes one or more clusters and
/// the clusters are appended in block order, so a coordinate-sorted input gives
/// a file that region queries can seek in. \p block_bytes is the size of the
/// input blocks handed to the workers (the default is a good value; tests use a
/// small one to exercise many blocks). The output holds the same records as
/// samtoramntuple() would write.
void samtoramntuple_parallel(const char *datafile, const char *treefile, int compression_algorithm,
                             uint32_t quality_policy, int threads, size_t block_bytes = 64u << 20);

/// Writes one RAM file per reference, named <output_prefix>_<rname>.root.
/// Records stream to their file as they are parsed, in the order they arrive.
void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy);

#endif
