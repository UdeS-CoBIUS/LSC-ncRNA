cmake_minimum_required(VERSION 3.15)
project(MotifsExtractionSelection)

set(CMAKE_CXX_STANDARD 17)

add_executable(MotifsExtractionSelection
        Main.cpp
        Node.h
        SuffixTree_QuadraticTime.cpp SuffixTree_QuadraticTime.h
        PrintTree.h
        Old_CommonMotifs.cpp Old_CommonMotifs.h
        FastaFilesReader.cpp FastaFilesReader.h
        SubSetDistancePercentage.h
        CommonMotifs.cpp CommonMotifs.h hasher_Container.h)