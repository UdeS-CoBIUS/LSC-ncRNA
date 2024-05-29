cmake_minimum_required(VERSION 3.15)
project(constructTrainTestFiles)

set(CMAKE_CXX_STANDARD 17)

add_executable(constructTrainTestFiles
        main.cpp
        FastaFilesReader.h FastaFilesReader.cpp
        Utile.h
        )