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

ROOT_EXECUTABLE(bamtoramntuple
    bamtoramntuple.cxx
    LIBRARIES
        ramcore
        hts
        ROOT::Core
        ROOT::ROOTNTuple
)

install(TARGETS samtoramntuple ramntupleview bamtoramntuple
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
