How the code works
==================

One shared library, ``ramcore``, holds everything; the command-line tools
are thin wrappers over it. This page follows a record from SAM text into
the file, a query from a region string to the rows it returns, and then
the one piece of state that everything shares.

Layout
------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Directory
     - Contents
   * - ``inc/rntuple``
     - ``RAMNTupleRecord``: the record, its encoders, the name tables, the
       index, and the metadata reader and writer
   * - ``inc/ramcore``
     - the SAM parser, the converters, the region scan
   * - ``src``
     - their implementations, same subdirectories
   * - ``tools``
     - the four command-line tools, one file each
   * - ``test``
     - googletest suites and ``ctest`` entries
   * - ``benchmark``
     - Google Benchmark binaries and the SAM generator
   * - ``scripts``
     - benchmark runner and table renderer
   * - ``docs``
     - this site; the API pages are generated from the headers

How a record gets in
--------------------

.. code-block:: text

   SAM line  ->  SamParser  ->  SamRecord  ->  RAMNTupleRecord setters  ->  RNTuple entry  ->  Fill()

1. **Read.** ``SamParser::ParseFile`` reads the input with ``getline``, so a
   line of any length works. Header lines go to one callback as
   ``(tag, content)``.
2. **Validate.** Each alignment line is split on tabs without collapsing
   empty fields. Every mandatory column is checked: integers in range, a
   well-formed CIGAR, nothing empty. A bad line is reported with its number
   and skipped; it never reaches the writer.
3. **Encode.** The converter copies the ``SamRecord`` into a
   ``RAMNTupleRecord`` through its setters. That is where the encoding
   happens: names become table indices, positions become 0-based, the CIGAR
   and the sequence are packed, and quality is stored according to the
   policy flag.
4. **Write.** The converter fills the RNTuple entry and, for a mapped
   record, decides whether this row gets an index entry (first on its
   reference, 10 kb past the last entry, or every 100th mapped read).
5. **Finish.** After the last record it writes ``INDEX``, ``METADATA`` and
   the ``headers`` key.

``bamtoramntuple`` follows the same steps from htslib's ``bam1_t``,
formatting each field to the text the SAM path would have seen, so both
produce identical files.

How a query gets out
--------------------

.. code-block:: text

   "chr1:1000-2000"  ->  ramntuplescan  ->  callback per matching row

``ramntuplescan`` in ``RAMNTupleView.cxx`` is the one region walk. Both
tools are built on it: ``ramntupleview`` passes a counting callback,
``ramdump`` passes one that reads the full record and formats it. If you
change what "in the region" means, change it there and both tools follow.

1. **Resolve the reference.** The whole string is tried as a name first,
   then split on the last colon, the order samtools uses.
2. **Seek.** On a sorted file the index gives the row at or before
   ``start - max_ref_span``. The backoff is what makes the sparse index
   exact: no read that begins earlier and reaches into the region can be
   skipped. On an unsorted file there is no seek.
3. **Scan.** Rows are read forward through the ``refid``, ``pos`` and
   ``cigar`` columns only. On a sorted file the scan stops at the first row
   past the region; on an unsorted file it runs to the end.
4. **Test.** A row is in the region when its reference matches and
   ``[pos, pos + span - 1]`` meets the query interval, htslib's
   ``bam_endpos`` rule. Secondary and supplementary reads count, as in
   ``samtools view``.
5. **Deliver.** Each matching row number goes to the callback; the caller
   owns whatever view it reads the record through.

Shared state
------------

Four things a query needs are not in the records. All four are static
members of ``RAMNTupleRecord``, one copy per process.

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - State
     - Filled by the writers
     - Loaded by ``OpenRAMFile()``
   * - RNAME and RNEXT name tables
     - from ``@SQ`` lines and records
     - from ``METADATA``
   * - region index
     - one entry per index rule hit
     - from ``INDEX``
   * - longest reference span
     - running maximum over records
     - from ``METADATA``
   * - sort flag
     - running order check
     - from ``METADATA``

``InitializeRefs()`` is the "start a new file" step: it creates the tables
if they do not exist and wipes the per-file parts, so a second file in the
same process does not inherit the first one's index, span or sort flag.
It has exactly four callers: ``samtoramntuple``, ``bamtoramntuple``,
``samtoramntuple_split_by_chromosome`` (once, before the first record),
and ``OpenRAMFile()``.

.. important::

   Constructing a record must never touch the shared state.

   RNTuple constructs a ``RAMNTupleRecord`` whenever a view of the
   ``record`` field is created and whenever a writer model is created, and
   both happen after the file state is loaded. The constructor calls only
   ``EnsureTables()``, which creates the tables and changes nothing.

What went wrong when the rule was broken
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The constructor once called ``InitializeRefs()``. Two symptoms followed:

- ``ramdump`` creates a record view after opening the file. The reset put
  the sort flag back to "sorted", the scan stopped at the first row past
  the region, and on an unsorted file records were dropped. The count-only
  path, which creates no view, was right, so the two disagreed.
- The split writer opens a model for each chromosome. The reset zeroed the
  longest span every time, so earlier chromosomes' files were written with
  a span too small for the seek to back off enough. A span that is too
  large costs time; one that is too small loses reads.

When you add per-file state
~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Put its reset in ``InitializeRefs()`` and nowhere else.
- Extend ``ConstructingARecordKeepsTheOpenFileState`` in
  ``test/ramcoretests.cxx``. It opens a file, creates a record view and a
  record, and checks the flag, the span and the index size are unchanged.
- If a writer computes it, extend ``SplitFilesKeepTheLongestSpan`` in
  ``test/chromosome_split_test.cxx`` the same way.
- Test ``ramdump`` on unsorted input, not only ``ramntupleview``; only
  ``ramdump`` goes through a record view.
- After writing split files, read the metadata back. Byte-identical dumps
  prove the records, not the metadata.
