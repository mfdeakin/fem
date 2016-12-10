
cmake_minimum_required (VERSION 3.3)
project(numerical)

include_directories(./include)

add_executable(tester src/test/test.cpp)

set_target_properties(tester PROPERTIES COMPILE_FLAGS "-g -std=c++14")
