Contributing
============

Layout
------

==================  ==========================================================
``inc/rntuple``     ``RAMNTupleRecord``: the record, its encoders, the name
                    tables, the index, and the metadata reader and writer
``inc/ramcore``     the SAM parser, the converters, the region scan
``src``             their implementations
``tools``           the four command-line tools, one file each
``test``            googletest suites and ``ctest`` entries
``benchmark``       Google Benchmark binaries and the SAM generator
``scripts``         benchmark runner and table renderer
``docs``            this site; the API pages are generated from the headers
==================  ==========================================================

How a record gets in
--------------------

``SamParser::ParseFile`` reads the input line by line with ``getline``, so
lines of any length work. Header lines go to one callback as ``(tag,
content)``. Each alignment line is split on tabs without collapsing empty
fields, every mandatory column is validated, and a valid record goes to the
second callback as a ``SamRecord``. Invalid records are reported and never
reach the writer.

The converter copies the ``SamRecord`` into a ``RAMNTupleRecord`` through
its setters, which do the encoding: names become table indices, positions
become 0-based, the CIGAR is packed, the sequence is packed, and quality is
encoded according to the policy flag. It then fills the RNTuple entry and,
for mapped records, decides whether this row gets an index entry. After the
last record it writes ``INDEX``, ``METADATA`` and ``headers``.

``bamtoramntuple`` does the same from htslib's ``bam1_t``, formatting each
field to the text the SAM path would have seen so both produce identical
files.

How a query gets out
--------------------

``ramntuplescan`` in ``RAMNTupleView.cxx`` is the one region walk. It
resolves the reference name, seeks with the index, and calls a callback
with each matching row. ``ramntupleview`` passes a counting callback;
``ramdump`` passes one that reads the full record and formats it. If you
change what "in the region" means, change it there and both tools follow.

Shared state, and the one rule about it
---------------------------------------

Four things a query needs are not in the records: the RNAME and RNEXT name
tables, the region index, the longest reference span, and the sort flag.
All four are static members of ``RAMNTupleRecord``, one copy per process.
``OpenRAMFile()`` loads them from a file's ``METADATA`` and ``INDEX``; the
writers fill them while converting and write them out at the end.

``InitializeRefs()`` is the "start a new file" step: it creates the tables
if they do not exist and wipes the per-file parts, so a second file in the
same process does not inherit the first one's index, span or sort flag.
Exactly three callers are allowed: ``samtoramntuple``, ``bamtoramntuple``,
and ``OpenRAMFile()``. ``samtoramntuple_split_by_chromosome`` calls it once,
before the first record.

**The rule: constructing a record must never touch this state.** RNTuple
constructs a ``RAMNTupleRecord`` whenever a view of the ``record`` field is
created and whenever a writer model is created, and both happen after the
file state is loaded. The constructor once called ``InitializeRefs()``, and
the two symptoms were exactly that:

- ``ramdump`` creates a record view after opening the file. The reset put
  the sort flag back to "sorted", the scan stopped at the first record past
  the region, and on an unsorted file records were dropped. The count-only
  path, which creates no view, was right, so the two disagreed.
- The split writer opens a model for each chromosome. The reset zeroed the
  longest span every time, so earlier chromosomes' files were written with
  a span too small for the region seek to back off enough. A span that is
  too large costs time; one that is too small loses reads.

The constructor now calls only ``EnsureTables()``, which creates the tables
and changes nothing. If you add per-file state, put its reset in
``InitializeRefs()`` and nowhere else, and extend
``ConstructingARecordKeepsTheOpenFileState`` in ``ramcoretests.cxx``: it
opens a file, creates a record view and a record, and checks the flag, the
span and the index size are unchanged. ``SplitFilesKeepTheLongestSpan`` in
``chromosome_split_test.cxx`` guards the writer side.

Two habits would have caught this earlier. Test ``ramdump`` on unsorted
input, not only ``ramntupleview``, because only ``ramdump`` goes through a
record view. And after writing split files, read their ``max_ref_span``
back; byte-identical dumps prove the records, not the metadata.

Tests
-----

``test/ramcoretests.cxx`` covers the parser, the encoders, the index rules
and region queries with small hand-written SAM inputs whose expected counts
are stated in the test. ``test/ramdump_test.cxx`` runs the built tools on a
fixture and checks the dump reproduces it and that counts follow samtools'
rules. ``test/chromosome_split_test.cxx`` checks split output against the
unsplit file. When you fix a bug, add the input that exposed it.

Run everything with ``ctest`` from the build directory. A single suite runs
faster on its own, for example ``build/test/ramcoretests
--gtest_filter='*Cigar*'``.

Style and CI
------------

- Code is formatted with clang-format 18 using the repository's
  ``.clang-format``; the ``precheckin`` job runs ``git clang-format``
  against the PR base and fails on any difference. Format only the lines
  you touch: ``git clang-format develop``.
- ``.clang-tidy`` is strict. The review job posts its findings as PR
  comments. New code should be clean; findings on lines you did not change
  are not yours to fix in the same PR.
- ``Build and Test`` runs ``ctest`` on Ubuntu with the ROOT release named in
  the CI workflow. ``Code Coverage`` posts to Codecov and fails a PR whose
  patch lowers coverage.

Pull requests
-------------

Branch from ``develop`` and open the PR against it. One topic per PR; one
change per commit with a message that says what was wrong and what changed,
in plain sentences. Put the evidence in the PR description: what you ran,
what it printed, what you compared it against. Reviewers here check claims
against samtools, so give them the command that shows the agreement.
