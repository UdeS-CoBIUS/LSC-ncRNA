cmake_minimum_required(VERSION 3.15)
project(testInfernal)

set(CMAKE_CXX_STANDARD 14)

add_executable(testInfernal testInfernal.cpp testInfernal_hmmonly.cpp)