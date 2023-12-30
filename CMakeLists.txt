cmake_minimum_required(VERSION 3.10)
project(GeneticAlgorithms)
add_executable(ga main.cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

target_include_directories(ga PRIVATE genetic_algorithm)
