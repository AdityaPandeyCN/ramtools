#ifndef RAMCORE_SAMTONTUPLE_H
#define RAMCORE_SAMTONTUPLE_H

#include <cstdint>
#include <string>

/// Converts a SAM file to one RAM file.
///
/// - `datafile`: input SAM
/// - `treefile`: output RAM file, created or overwritten
/// - `index`: write the sparse region index
/// - `split`, `cache`: accepted but not used; splitting is samtoramntuple_split_by_chromosome()
/// - `compression_algorithm`: ROOT compression code, algorithm * 100 + level (505 is ZSTD level 5)
/// - `quality_policy`: one of RAMNTupleRecord::EQualCompressionBits
///
/// Malformed records are reported on standard error and skipped. Prints the
/// reference tables and the record and index counts when done.
void samtoramntuple(const char *datafile,
                    const char *treefile,
                    bool index, bool split, bool cache,
                    int compression_algorithm,
                    uint32_t quality_policy);

/// Converts a SAM file to one RAM file per reference sequence, named
/// `<output_prefix>_<rname>.root`. Records with no reference are not written.
/// \a num_threads bounds the writer threads used per file.
void samtoramntuple_split_by_chromosome(const char *datafile, const char *output_prefix, int compression_algorithm,
                                        uint32_t quality_policy, int num_threads = 4);

#endif
