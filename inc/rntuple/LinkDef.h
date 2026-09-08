// NOLINTBEGIN(llvm-header-guard)
// ROOT dictionary for rootcling. Not a regular header: an include guard wraps
// the pragmas and can break dictionary generation.
#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RAMNTupleRecord + ;
#pragma link C++ class RAMNTupleIndex + ;
// IndexEntry is stored in every RAM file, so ROOT needs a TClass for it on each
// open. Without a dictionary entry cling parses the headers instead, which costs
// about two seconds per open.
// clang-format off
#pragma link C++ class RAMNTupleIndex::IndexEntry + ;
#pragma link C++ class std::vector<RAMNTupleIndex::IndexEntry> + ;
// clang-format on

#endif
// NOLINTEND(llvm-header-guard)
