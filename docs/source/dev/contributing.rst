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
