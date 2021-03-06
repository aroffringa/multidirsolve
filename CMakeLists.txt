cmake_minimum_required(VERSION 2.8.8)

project(mds)

add_definitions(-DAOPROJECT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -ggdb -Wall -march=native -std=c++11 -fopenmp")

#if(NOT DISABLE_OMP)
#  FIND_PACKAGE( OpenMP REQUIRED)
#  if(OPENMP_FOUND)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#  endif(OPENMP_FOUND)
#endif(NOT DISABLE_OMP)

# Casacore has a separate CMake file in this directory
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake)

SET(CASACORE_MAKE_REQUIRED_EXTERNALS_OPTIONAL TRUE)

# The following stuff will set the "rpath" correctly, so that
# LD_LIBRARY_PATH doesn't have to be set.

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

find_package(Threads REQUIRED)
#find_library(ARMADILLO_LIB NAMES armadillo)
#find_path(ARMADILLO_INCLUDE_DIR NAMES armadillo)
#include_directories(${ARMADILLO_INCLUDE_DIR})

add_executable(mds mds.cc MultiDirSolver.cc Constraint.cc PhaseFitter.cc SmoothnessConstraint.cc TECConstraint.cc Stopwatch.cc)
target_link_libraries(mds ${CMAKE_THREAD_LIBS_INIT} lapack)

add_executable(fittec fittec.cc PhaseFitter.cc)

add_executable(kerneltest kerneltest.cc)
