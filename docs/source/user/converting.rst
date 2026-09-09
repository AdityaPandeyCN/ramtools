Converting to RAM
=================

Two tools write RAM files: ``samtoramntuple`` reads SAM text and
``bamtoramntuple`` reads BAM through htslib. They take the same options and
produce the same file layout.

.. code-block:: bash

   samtoramntuple <input.sam> [output] [options]
   bamtoramntuple <input.bam> [output] [options]

If ``output`` is omitted the input name is used with its extension removed.
Unless the name already ends in ``.ram`` or ``.root``, ``.ram`` is appended.
A RAM file is an ordinary ROOT file; either extension works.

Options
-------

``-noindex``
   Do not write the region index. Region queries on the file then scan from
   the first record of the reference. Use it for files you will only ever
   read end to end.

``-illumina``
   Store quality scores with Illumina's 8-level binning instead of the full
   Phred range. Every score maps to one of 0, 1, 6, 15, 22, 27, 33, 37 or 40,
   which compresses far better and is what most variant callers are
   evaluated against. The mapping is lossy.

``-dropqual``
   Store no quality scores at all. Reads come back with ``*`` in the QUAL
   column.

``-compression N``
   The ROOT compression code, ``algorithm * 100 + level``. The default is
   505, ZSTD level 5.

   ========  =========  ======================================================
   Code      Algorithm  When to use it
   ========  =========  ======================================================
   ``505``   ZSTD 5     the default: good ratio, fast decompression
   ``404``   LZ4 4      fastest reads, larger files
   ``101``   ZLIB 1     compatibility with tools that lack ZSTD
   ``207``   LZMA 7     smallest files, slow to write and to read
   ``0``     none       for measuring the cost of compression itself
   ========  =========  ======================================================

   Algorithms 1, 2, 4 and 5 with levels 1 to 9 are accepted; anything else is
   rejected with a message. ``samtoramntuple`` only:

``-split``
   Write one file per reference sequence instead of one file, named
   ``<output>_<rname>.root``. Records with no reference (``*``) are not
   written. See :ref:`split`.

What the input has to look like
-------------------------------

**Sorted input gets an index.** The index lets a region query jump close to
the region and stop at the first record past it. Both steps assume the
records are in coordinate order, the order ``samtools sort`` produces. If
your input is not sorted, sort it first, or the index will send queries to
the wrong place:

.. code-block:: bash

   samtools sort -o sorted.bam reads.bam
   bamtoramntuple sorted.bam reads.ram

**Malformed records are skipped, not stored.** The SAM parser validates every
mandatory column: integer fields must be integers in range, the CIGAR must be
``*`` or a run of ``<length><op>`` pairs with the operators ``MIDNSHP=X``,
and no mandatory column may be empty. A record that fails is reported with
its line number on standard error and left out. The counts printed at the
end tell you how many records were written.

**Everything else round-trips.** Read names, FLAG, positions, MAPQ, CIGAR,
mate fields, TLEN, sequence, quality and every optional tag are stored so
that ``ramdump -h`` reproduces the input line for line. Lower-case bases are
stored as their upper-case code; a base outside ``ACGTN`` and the IUPAC
codes is stored as ``N``.

.. _split:

Splitting by chromosome
-----------------------

``samtoramntuple reads.sam out -split`` writes ``out_chr1.root``,
``out_chr2.root`` and so on, one complete RAM file per reference sequence,
each with its own metadata and a copy of the header. Query and dump them
like any other RAM file. Records with no reference sequence are not written
to any of them.

What gets printed
-----------------

Both tools print the reference tables, the number of records written and,
with an index, the number of index entries. ``samtoramntuple`` also prints
how many header lines and records it read, so the difference from the
number written is the number of skipped records.
