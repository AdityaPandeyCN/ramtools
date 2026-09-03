#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstdint>
#include <string>

void samtoramntuple(const char *datafile,
                    const char *treefile,
                    bool index, bool split, bool cache,
                    int compression_algorithm,
                    uint32_t quality_policy);

/// Writes one RAM file per reference, named `<output_prefix>_<rname>.root`.
///
/// Records stream to their chromosome's file as they are parsed, so memory is a
/// few page buffers per open file rather than the whole input, and each file
/// keeps the order its records arrived in and gets its own index.
void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy);

#endif
