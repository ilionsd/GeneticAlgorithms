cmake_minimum_required(VERSION 3.10)
project(GeneticAlgorithms)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

file(GLOB SOURCES
    genetic_algorithm/*.h
    main.cpp
)

add_executable(GeneticAlgorithms ${SOURCES})
