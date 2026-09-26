ROOT_EXECUTABLE(samtoramntuple
    samtoramntuple.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(ramntupleview
    ramntupleview.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(ramdump
    ramdump.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::RIO
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(bamtoramntuple
    bamtoramntuple.cxx
    LIBRARIES
        ramcore
        hts
        ROOT::Core
        ROOT::ROOTNTuple
)

ROOT_EXECUTABLE(raminfo
    raminfo.cxx
    LIBRARIES
        ramcore
        ROOT::Core
        ROOT::RIO
        ROOT::ROOTNTuple
        ROOT::ROOTNTupleUtil
)

install(TARGETS samtoramntuple ramntupleview ramdump bamtoramntuple raminfo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
