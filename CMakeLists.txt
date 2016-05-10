CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(groupAnalysis)

FIND_PACKAGE(statismo REQUIRED)
FIND_PACKAGE(VTK REQUIRED)
include (${VTK_USE_FILE})

include_directories(${statismo_INCLUDE_DIRS})
link_directories( ${VTK_LIBRARY_DIRS} ${statismo_LIBRARY_DIRS})

add_executable (groupAnalysis groupAnalysis.cxx)
target_link_libraries (groupAnalysis ${VTK_LIBRARIES} ${statismo_LIBRARIES})
