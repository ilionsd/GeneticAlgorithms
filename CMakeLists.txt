cmake_minimum_required(VERSION 3.10)
project(GeneticAlgorithms)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)

file(GLOB SOURCES
#    io/*.h
    genetic_algorithm/common/*.h
    genetic_algorithm/*.h
    *.h
    *.cpp
)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(GeneticAlgorithms ${SOURCES})
target_include_directories(GeneticAlgorithms PUBLIC ${CMAKE_SOURCE_DIR}/third-party)
target_link_libraries(GeneticAlgorithms ${Boost_LIBRARIES})
