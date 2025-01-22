# include ROOT
include_directories(${ROOT_INCLUDE_DIRS})

ROOT_GENERATE_DICTIONARY(G__GEScaleSyst
    GEScaleSyst.h
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/GEScaleSyst_LinkDef.h
)

add_library(GEScaleSyst SHARED
    GEScaleSyst.cc
    G__GEScaleSyst.cxx
)

target_include_directories(GEScaleSyst PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(GEScaleSyst ${ROOT_LIBRARIES})

install(TARGETS GEScaleSyst DESTINATION lib)

set(ROOTMAP "${CMAKE_CURRENT_BINARY_DIR}/libGEScaleSyst.rootmap")
set(PCM "${CMAKE_CURRENT_BINARY_DIR}/libGEScaleSyst_rdict.pcm")
message(STATUS "ROOTMAP: ${ROOTMAP}")
message(STATUS "PCM: ${PCM}")
install(FILES ${ROOTMAP} ${PCM} DESTINATION lib)
