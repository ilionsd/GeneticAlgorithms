cmake_minimum_required(VERSION 3.10)
project(GeneticAlgorithms)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB SOURCES
    genetic_algorithm/common/*.h
    genetic_algorithm/*.h
    *.h
    *.cpp
)

add_executable(GeneticAlgorithms ${SOURCES})

target_link_libraries(GeneticAlgorithms ${Boost_LIBRARIES})

target_include_directories(GeneticAlgorithms PUBLIC ${CMAKE_SOURCE_DIR}/third-party/jsoncpp/include)
target_link_libraries(GeneticAlgorithms ${CMAKE_SOURCE_DIR}/third-party/jsoncpp/lib/libjsoncpp.a)
