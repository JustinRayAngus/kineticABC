project("TDBoltzmann")
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# include stuff for HDF5
set(ENV{HDF5_HOME} "/Users/angus1/HDF5/hdf5-1.10.0-patch1/hdf5")
#find_package(HDF5 "1.8.14" REQUIRED COMPONENTS C HL)
#find_package(HDF5 REQUIRED COMPONENTS C HL)
#include_directories(${HDF5_INCLUDE_DIRS})
#link_directories(${HDF5_LIBRARY_DIRS})
include_directories($ENV{HDF5_HOME}/include)
link_directories($ENV{HDF5_HOME}/lib)

add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

set(sources 
   main.cpp
   jsoncpp.cpp)
add_executable(tdb ${sources})
target_link_libraries(tdb hdf5 hdf5_cpp)
