Name tables and index
=====================

Both classes are owned by ``RAMNTupleRecord`` as static members and are
serialised into the ``METADATA`` and ``INDEX`` ntuples. After
``OpenRAMFile()`` you reach them through ``GetRnameRefs()``,
``GetRnextRefs()`` and ``GetIndex()``.

RAMNTupleRefs
-------------

.. doxygenclass:: RAMNTupleRefs

.. doxygenfunction:: RAMNTupleRefs::GetRefId

.. doxygenfunction:: RAMNTupleRefs::FindRefId

.. doxygenfunction:: RAMNTupleRefs::GetRefName

.. doxygenfunction:: RAMNTupleRefs::Size

.. doxygenfunction:: RAMNTupleRefs::GetRefs

.. doxygenfunction:: RAMNTupleRefs::AddRef

.. doxygenfunction:: RAMNTupleRefs::SetRefs

.. doxygenfunction:: RAMNTupleRefs::Clear

RAMNTupleIndex
--------------

.. doxygenclass:: RAMNTupleIndex

.. doxygenstruct:: RAMNTupleIndex::IndexEntry
   :members:

.. doxygenfunction:: RAMNTupleIndex::AddItem

.. doxygenfunction:: RAMNTupleIndex::GetRow

.. doxygenfunction:: RAMNTupleIndex::GetRowsInRange

.. doxygenfunction:: RAMNTupleIndex::Size

.. doxygenfunction:: RAMNTupleIndex::GetEntries

.. doxygenfunction:: RAMNTupleIndex::SetEntries

.. doxygenfunction:: RAMNTupleIndex::Clear
