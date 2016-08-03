project("TDBoltzmann")
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# include stuff for HDF5
find_package(HDF5 "1.8.14" REQUIRED COMPONENTS C HL)
include_directories(${HDF5_INCLUDE_DIRS})
link_directories(${HDF5_LIBRARY_DIRS})

add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set(sources 
   main.cpp
   jsoncpp.cpp)
add_executable(tdb ${sources})
target_link_libraries(tdb hdf5 hdf5_cpp)
