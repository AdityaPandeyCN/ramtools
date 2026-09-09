Benchmarks
==========

Four Google Benchmark binaries are built with ``RAMTOOLS_BUILD_BENCHMARKS``
and live in ``build/benchmark``:

======================================  ======================================
``sam_to_ram_benchmark``                conversion throughput on generated SAM
``conversion_time_benchmark``           conversion time by input size
``chromosome_split_benchmark``          ``-split`` against ``samtools`` splitting
``region_query_benchmark``              region queries on a file you provide
======================================  ======================================

The first three generate their own input. The region query benchmark
measures a real file, so it takes it from the environment:

.. code-block:: bash

   export RAMTOOLS_BENCH_RNTUPLE=/data/HG00154.ram
   export RAMTOOLS_BENCH_REGIONS="chr1:1000000-1001000,chr1:1-50000000,chr21"
   build/benchmark/region_query_benchmark

``RAMTOOLS_BENCH_REGIONS`` is optional; without it a fixed set of ``chr1``
and ``chr21`` regions is used. Each region is one benchmark case, and the
label of each case shows the region and how many records it returned, so a
zero is easy to spot. The binary refuses to run without
``RAMTOOLS_BENCH_RNTUPLE`` rather than measure an open that failed.

``scripts/run_benchmark.py`` runs the set and ``scripts/render_benchmark.py``
turns the JSON output into the tables in the README.

Reproducing a number
--------------------

A benchmark result is only a result if someone else can get it. Record
alongside every number: the input file and how it was made (``samtools
sort`` order, index or not), the compression code, the ROOT version, and
the machine. The README's tables were measured on the HG00154 sample from
the 1000 Genomes Project, 196 million reads, converted at LZMA, LZ4 and
ZLIB.
