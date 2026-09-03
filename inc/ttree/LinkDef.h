#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class RAMRefs+;
#pragma link C++ class RAMIndex+;
#pragma link C++ class RAMRecord+;
#pragma link C++ class RAMNTupleRecord+;
#pragma link C++ class RAMNTupleIndex+;
// IndexEntry is written into every RAM file, so ROOT builds a TClass for it the
// moment one is opened. Without these two lines there is no dictionary for it
// and the fallback is cling parsing the headers -- two seconds on every single
// open, dwarfing the query it was opened for.
#pragma link C++ class RAMNTupleIndex::IndexEntry+;
#pragma link C++ class std::vector<RAMNTupleIndex::IndexEntry>+;

#endif

