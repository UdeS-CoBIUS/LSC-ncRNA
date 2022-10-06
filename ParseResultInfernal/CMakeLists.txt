cmake_minimum_required(VERSION 3.15)
project(ParseResultInfernal)

set(CMAKE_CXX_STANDARD 17)

add_executable(ParseResultInfernal
        main.cpp
        FastaFilesReader.h
        FastaFilesReader.cpp
        Utile.h
        )