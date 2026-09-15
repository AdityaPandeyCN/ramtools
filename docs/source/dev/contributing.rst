Contributing
============

How the code is put together is on :doc:`architecture`; this page is about
getting a change in.

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
