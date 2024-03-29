
# CMakeLists.txt for FAS package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(MPDFlowAnalysis)

# You need to tell CMake where to find the ROOT installation. This can be done
# in a number of ways:
#   - ROOT built with classic configure/make use the provided
#   $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix
#   for ROOT

find_package(Git)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net MathMore EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

add_definitions(${ROOT_CXX_FLAGS})

add_subdirectory(PicoDst)
add_subdirectory(Readers)
add_subdirectory(StFemtoDst)
add_subdirectory(FlowAnalysis)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")
set (CMAKE_CXX_STANDARD 11)

set(FLOW_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/PicoDst/include
  ${CMAKE_CURRENT_SOURCE_DIR}/Readers
  ${CMAKE_CURRENT_SOURCE_DIR}/StFemtoDst
  ${CMAKE_CURRENT_SOURCE_DIR}/FlowAnalysis
  ${ROOT_INLCUDE_DIRS}
)

include_directories(${FLOW_INCLUDE_DIRECTORIES})

set(FLOW_INCLUDE_LIBRARIES
  StFemtoDst
  PicoDst
  Readers
  FlowAnalysis
  ${ROOT_LIBRARIES}
)

# set_source_files_properties(FlowQCumulant.C PROPERTIES LANGUAGE CXX)

add_executable(RunFlowAnalysis RunFlowAnalysis.cpp)
target_link_libraries(RunFlowAnalysis ${FLOW_INCLUDE_LIBRARIES})
add_executable(FemtoDstAnalyzer FemtoDstAnalyzer.cpp)
target_link_libraries(FemtoDstAnalyzer ${FLOW_INCLUDE_LIBRARIES})
add_executable(GetDCA GetDCA.cpp)
target_link_libraries(GetDCA ${FLOW_INCLUDE_LIBRARIES})
add_executable(FitDCA FitDCA.cpp)
target_link_libraries(FitDCA ${FLOW_INCLUDE_LIBRARIES})

add_executable(Recentering RunRecentering.cpp)
target_link_libraries(Recentering ${FLOW_INCLUDE_LIBRARIES})
add_executable(Flattening RunFlattening.cpp)
target_link_libraries(Flattening ${FLOW_INCLUDE_LIBRARIES})
