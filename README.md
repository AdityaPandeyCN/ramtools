# RAMTools - ROOT Alignment/Map Format Tools

RAMTools provides efficient tools for converting SAM files to ROOT's modern, columnar RNTuple format (RAM - ROOT Alignment/Map) and working with genomic alignment data.

## Features

- High-performance SAM to RAM conversion
- Chromosome-based splitting for parallel processing
- Region-based querying capabilities

## Requirements

- ROOT 6.38+
- C++17 compatible compiler
- CMake 3.16+

## Quick Start
```bash
# 1. Build the tools
mkdir build && cd build
cmake ..
make -j$(nproc)

# 2. Convert a SAM file to the RAM format
./tools/samtoramntuple input.sam output.root

# 3. Query a specific region from the command line
./tools/ramntupleview output.root "chr1:15700-15800"
```

## Command-Line Tools

The primary way to interact with RAMTools is through these command-line executables.

### SAM to RAM Conversion

Convert a standard SAM file into the optimized RNTuple-based RAM format.
```bash
# Basic conversion
./tools/samtoramntuple input.sam output.root

# Split by chromosome for parallel processing
# (Creates output_chr1.root, output_chr2.root, etc.)
./tools/samtoramntuple input.sam output -split
```

Options: `-noindex` skips the region index, `-illumina` stores 8-level binned
quality scores, `-dropqual` stores none, `-compression N` sets the ROOT
compression code (algorithm*100+level; the default 505 is ZSTD level 5).

The index needs the input in coordinate order. An unsorted input converts
fine but gets no index, and region queries on it read the whole file.

### Region Querying

Count the records overlapping a genomic region, as `samtools view -c` does.
```bash
# Usage: ./tools/ramntupleview [input.root] "[chromosome]:[start]-[end]"
./tools/ramntupleview output.root "chr1:10150-10300"
```

### Dumping back to SAM

`ramdump` writes a RAM file out as SAM and takes the `samtools view` options
`-h`, `-H`, `-c`, `-f`, `-F` and `-o`, so its output can be checked against
samtools directly.
```bash
# The whole file with its header; should reproduce the input SAM
./tools/ramdump -h output.root > roundtrip.sam

# Count primary alignments in a region
./tools/ramdump -c -F 0x900 output.root "chr1:10150-10300"
```

## Benchmarks

HG00154 from the 1000 Genomes Project: 196,040,370 records, 72.06 GB as SAM, coordinate sorted, GRCh37 reference names (`1`, not `chr1`). Run on 4 cores with 11 GB of memory, input and output on the same hard disk, samtools 1.13 with 4 threads, `samtoramntuple` on one thread. Every output holds the same 196,040,370 records. Query times are for the binary-search seek of #83.

### File size

![File size of HG00154 in each format](assets/benchmark_sizes.svg)

BAM is `samtools view -b`, CRAM is `samtools view -C -T ref.fasta`, and each RAM file is `samtoramntuple -compression N` with the flag shown. A CRAM cannot be read without its reference, so the stacked bar is the size of a self-contained copy. Uncompressed RAM (`-compression 0`) is 83.7 GB, larger than the SAM, because every string and vector column carries an offset column that the codecs then remove. ZSTD 9 is the only setting that comes in under CRAM plus its reference, at ten times the conversion time of ZSTD 5 for 9% less space.

### Region queries

Records overlapping a region, counted with `samtools view -c` for BAM and CRAM and `ramdump -c` for RAM (ZSTD 5). Best of three runs, warm cache. All three formats return the same records for every region.

| Region | Records | BAM (s) | CRAM (s) | RAM (s) |
|--------|---------|---------|----------|---------|
| 100 bp, `1:1000000-1000100` | 2 | 0.15 | 0.10 | 0.48 |
| BRCA2, `13:32889611-32973805` | 6,057 | 0.10 | 0.09 | 0.39 |
| 10 Mb, `1:10000000-20000000` | 582,985 | 0.80 | 1.75 | 0.40 |
| 100 Mb, `2:1-100000000` | 7,048,385 | 8.59 | 21.43 | 0.85 |

A RAM query spends about 0.4 s starting ROOT and opening the file, then finds the region by binary search and counts about 15 million records per second, since a count reads only the position and CIGAR columns. It is slower than samtools on small regions and faster from about 10 Mb up.

### Compression flags

`-compression N` takes a ROOT compression code, algorithm times 100 plus level: 1 is ZLIB, 2 LZMA, 4 LZ4, 5 ZSTD, with levels 1 to 9, and 0 means no compression. The default is 505, ZSTD level 5. LZ4 converts fastest and gives the largest files, LZMA and high ZSTD levels the smallest at a much higher conversion cost; the codec makes little difference to query time.
