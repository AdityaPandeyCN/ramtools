# RAMTools

RAM (ROOT Alignment/Map) is a file format for aligned sequencing reads, the same
data as a BAM, stored column by column in ROOT's RNTuple format. A RAM file is
about a quarter smaller than the BAM of the same reads, needs no reference to be
read, and answers region queries such as `chr13:32889611-32973805` the way
`samtools view` does. RAMTools converts SAM and BAM files to RAM, queries them,
and writes them back out as SAM so every result can be checked against samtools.

The tools:

| Tool | What it does |
|------|--------------|
| `samtoramntuple` | converts a SAM file to RAM |
| `bamtoramntuple` | converts a BAM file to RAM |
| `ramdump` | writes a RAM file, or a region of it, out as SAM; takes the `samtools view` options |
| `ramntupleview` | counts the records in a region and reports the query time |

## Installation

RAMTools builds on Linux. The steps below are for Ubuntu; other distributions
need the same three things: ROOT, htslib and a C++ compiler with CMake.

### 1. System packages

```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git pkg-config libhts-dev libtbb-dev libvdt-dev
```

`libhts-dev` is htslib, which reads the BAM input. `libtbb-dev` and `libvdt-dev`
are needed by ROOT.

### 2. ROOT

RAMTools needs ROOT 6.38 or newer. The quickest way is the prebuilt archive from
the ROOT project; pick the one matching your Ubuntu version from
https://root.cern/install/all_releases/. For Ubuntu 24.04:

```bash
wget https://github.com/root-project/root/releases/download/v6-38-06/root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz
sudo tar -xzf root_v6.38.06.Linux-ubuntu24.04-x86_64-gcc13.3.tar.gz -C /opt/
```

ROOT has to be put on your path in every terminal you use it from:

```bash
source /opt/root/bin/thisroot.sh
```

Add that line to your `~/.bashrc` to make it permanent. `root --version` should
then print 6.38. If you use conda, `conda install -c conda-forge root` works as
well and needs no `thisroot.sh`.

### 3. Build RAMTools

```bash
git clone https://github.com/compiler-research/ramtools.git
cd ramtools
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

The tools are in `build/tools/`. Running one without arguments prints its
usage:

```bash
./build/tools/samtoramntuple
```

The tests are optional and take a minute:

```bash
ctest --test-dir build
```

## Converting to RAM

From a SAM file:

```bash
./build/tools/samtoramntuple reads.sam reads.ram
```

From a BAM file:

```bash
./build/tools/bamtoramntuple reads.bam reads.ram
```

Both print the reference names and the number of records written. Options for
either tool:

| Option | Effect |
|--------|--------|
| `-compression N` | ROOT compression code, algorithm times 100 plus level; the default 505 is ZSTD level 5 (see below) |
| `-illumina` | stores quality scores in Illumina's 8 bins, which makes the file smaller |
| `-dropqual` | stores no quality scores |
| `-threads N` | compresses on N threads; four threads roughly halve the conversion time |
| `-split` | `samtoramntuple` only: one RAM file per chromosome, `reads_chr1.root`, `reads_chr2.root`, and so on |

Sort the input by coordinate first, as `samtools sort` does. An unsorted file
converts fine, but region queries on it have to read the whole file instead of
jumping to the region; the converter says so when that happens.

## Querying a region

Regions are written as in samtools: `chr1`, `chr1:1000`, or `chr1:1000-2000`,
1-based and inclusive. Reference names are whatever the input used, so a
GRCh37 file is queried as `13:32889611-32973805`, not `chr13`.

Count the records overlapping a region, as `samtools view -c` does:

```bash
./build/tools/ramdump -c reads.ram chr1:10150-10300
```

Print them as SAM lines:

```bash
./build/tools/ramdump reads.ram chr1:10150-10300
```

`ramdump` takes the `samtools view` options `-h` (include the header), `-H`
(header only), `-c` (count only), `-f` and `-F` (keep records with all, or none,
of these FLAG bits) and `-o FILE`. For example, primary alignments only:

```bash
./build/tools/ramdump -c -F 0x900 reads.ram chr1:10150-10300
```

`ramntupleview` counts as well and prints how long the query took:

```bash
./build/tools/ramntupleview reads.ram chr1:10150-10300
```

## Checking a RAM file against samtools

A RAM file holds everything the input did, so it can be written back out and
compared:

```bash
./build/tools/ramdump -h reads.ram > roundtrip.sam
samtools view -c reads.bam chr1:10150-10300
./build/tools/ramdump -c reads.ram chr1:10150-10300
```

The two counts must agree, and `roundtrip.sam` reproduces the input line for
line.

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

## Contributing

Bug reports, test cases and pull requests are welcome. [CONTRIBUTING.md](CONTRIBUTING.md)
describes the development build, the test suites, the formatting and clang-tidy
checks that run on every pull request, and what a reviewer looks for. The short
version: branch from `develop`, keep each PR to one topic, and put the samtools
command that confirms your result in the description.

## Reporting a problem

Open an issue at https://github.com/compiler-research/ramtools/issues. Include
the command you ran, its output, your ROOT version, and for a wrong result the
samtools command that gives the answer you expected. A few SAM lines that
reproduce the problem are the most useful thing you can attach.

## License

RAMTools is released under the Apache License 2.0; see [LICENSE](LICENSE).
