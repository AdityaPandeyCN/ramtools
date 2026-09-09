Converters
==========

The writers. Each calls ``RAMNTupleRecord::InitializeRefs()``, fills records
through the setters, decides per record whether it gets an index entry, and
writes ``INDEX``, ``METADATA`` and the ``headers`` key when the input is
exhausted. The conversion options are described in :doc:`../user/converting`.

samtoramntuple
--------------

.. doxygenfunction:: samtoramntuple

samtoramntuple_split_by_chromosome
----------------------------------

.. doxygenfunction:: samtoramntuple_split_by_chromosome

bamtoramntuple
--------------

.. doxygenfunction:: bamtoramntuple
