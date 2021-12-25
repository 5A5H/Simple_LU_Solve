cmake_minimum_required(VERSION 3.0)
project(LUSOLVE VERSION 0.1)

set(CMAKE_CXX_STANDARD 17)

include(CTest)
enable_testing()

add_executable(lu_solve main.cpp)

add_test(lu_solve lu_solve)