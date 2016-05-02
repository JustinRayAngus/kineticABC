project(JustinsFirstCode)
cmake_minimum_required(VERSION 2.6)

include_directories($ENV{GMOCK_HOME}/include $ENV{GTEST_HOME}/include $ENV{JSONCPP_HOME}/include $ENV{CURL_HOME}/include)
link_directories($ENV{GMOCK_HOME}/build $ENV{GTEST_HOME}/build $ENV{JSONCPP_HOME}/build/src/lib_json $ENV{CURL_HOME}/build/lib)
add_definitions(-std=c++0x)
set(CMAKE_CXX_FLAGS "${CMAXE_CXX_FLAGS} -Wall")

set(sources 
   main.cpp)
add_executable(test ${sources})
target_link_libraries(test pthread)
target_link_libraries(test gmock)
target_link_libraries(test gtest)
target_link_libraries(test curl)
target_link_libraries(test jsoncpp)

