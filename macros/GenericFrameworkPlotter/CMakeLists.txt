cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PlotQCumulant)


# find_package(Git)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net EG MathMore)
include(${ROOT_USE_FILE})
add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")

set(ACCCOR_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOT_INLCUDE_DIRS}
)

include_directories(${ACCCOR_INCLUDE_DIRECTORIES})

set(FLOWANA_INCLUDE_LIBRARIES
  ${ROOT_LIBRARIES}
)

set(FlowANA_h_files
  "${CMAKE_CURRENT_SOURCE_DIR}/MultiparticleCorrelation.h"
)

set(FlowANA_cxx_files
  "${CMAKE_CURRENT_SOURCE_DIR}/MultiparticleCorrelation.cxx"
)

set(FlowANA_LinkDef
	"${CMAKE_CURRENT_SOURCE_DIR}/MultiparticleCorrelation.LinkDef.h"
)

#---Generate dictionary
ROOT_GENERATE_DICTIONARY(G__MultiparticleCorrelation
  ${FlowANA_h_files}
  LINKDEF ${FlowANA_LinkDef}
)

#---Compile library
add_library(MultiparticleCorrelation SHARED ${FlowANA_cxx_files} G__MultiparticleCorrelation.cxx)
target_link_libraries(MultiparticleCorrelation ${FLOWANA_INCLUDE_LIBRARIES})

set(ACCCOR_INCLUDE_LIBRARIES
  MultiparticleCorrelation
  ${ROOT_LIBRARIES}
)

add_executable(PlotV2QCumulant PlotV2QCumulant.C) # main.cpp
target_link_libraries(PlotV2QCumulant ${ACCCOR_INCLUDE_LIBRARIES})
