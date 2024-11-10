if (APPLE)
    set(LIB_EXT "dylib")
else()
    set(LIB_EXT "so")
endif()

# include ROOT and LHAPDF
set(LHAPDF_INCLUDE_DIR $ENV{LHAPDF_INCLUDE_DIR})
include_directories(${LHAPDF_INCLUDE_DIR})
include_directories(${ROOT_INCLUDE_DIRS})

# Generate the ROOT dictionary before defining the library
file(GLOB HEADERS include/*.h)
list(FILTER HEADERS EXCLUDE REGEX "AnalyzerTools_LinkDef.h")
ROOT_GENERATE_DICTIONARY(G__AnalyzerTools 
    ${HEADERS}
    LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/AnalyzerTools_LinkDef.h
)

# Define the library with source files and the generated dictionary
file(GLOB SOURCES src/*.C)
add_library(AnalyzerTools SHARED
    ${SOURCES}
    G__AnalyzerTools.cxx
)

# Specify where the library should look for header files
target_include_directories(AnalyzerTools PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link ROOT and lhapdf libraries
set(LHAPDF_LIBRARY $ENV{LHAPDF_LIB_DIR}/libLHAPDF.${LIB_EXT})
target_link_libraries(AnalyzerTools ${ROOT_LIBRARIES})
target_link_libraries(AnalyzerTools ${LHAPDF_LIBRARY})
target_link_libraries(AnalyzerTools DataFormats)

# Specify where to install the library
install(TARGETS AnalyzerTools DESTINATION lib)

# Iinstall .rootmap and .pcm files to the lib
set(ROOTMAP "${CMAKE_CURRENT_BINARY_DIR}/libAnalyzerTools.rootmap")
set(PCM "${CMAKE_CURRENT_BINARY_DIR}/libAnalyzerTools_rdict.pcm")
message(STATUS "ROOTMAP: ${ROOTMAP}")
message(STATUS "PCM: ${PCM}")
install(FILES ${ROOTMAP} ${PCM} DESTINATION lib)
