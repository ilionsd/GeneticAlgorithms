cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW)

project(GeneticAlgorithms)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

add_library(jsoncpp STATIC IMPORTED)
if(${CMAKE_SYSTEM_NAME} STREQUAL "Windows")
    set_target_properties(
        # Specifies the target library.
        jsoncpp

        # Specifies the parameter you want to define.
        PROPERTIES IMPORTED_LOCATION

        # Provides the path to the library you want to import.
        ${CMAKE_SOURCE_DIR}/third-party/jsoncpp/lib/jsoncpp.lib)
else()
    set_target_properties(
        # Specifies the target library.
        jsoncpp

        # Specifies the parameter you want to define.
        PROPERTIES IMPORTED_LOCATION

        # Provides the path to the library you want to import.
        ${CMAKE_SOURCE_DIR}/third-party/jsoncpp/lib/libjsoncpp.a)
endif()

set( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}
        "/NODEFAULTLIBS:LIBCMT" )
set( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}
        "/NODEFAULTLIBS:MSVCRT" )
set( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS}
        "/NODEFAULTLIBS:MSVCPRTD" )

file(GLOB SOURCES
        genetic_algorithm/dummy/*h
        genetic_algorithm/cmn/*.h
        genetic_algorithm/*.h
        *.h
        *.cpp
)

add_executable(GeneticAlgorithms ${SOURCES})

target_link_libraries(GeneticAlgorithms ${Boost_LIBRARIES})

target_include_directories(GeneticAlgorithms PUBLIC 
        ${CMAKE_SOURCE_DIR}/third-party/jsoncpp/include)

target_link_libraries(GeneticAlgorithms 
        jsoncpp
)

set_target_properties(GeneticAlgorithms PROPERTIES LINK_FLAGS 
        "/NODEFAULTLIB:user32.lib
        /NODEFAULTLIB:kernel32.lib
        /NODEFAULTLIB:gdi32.lib
        /NODEFAULTLIB:winspool.lib
        /NODEFAULTLIB:shell32.lib
        /NODEFAULTLIB:ole32.lib
        /NODEFAULTLIB:oleaut32.lib
        /NODEFAULTLIB:uuid.lib
        /NODEFAULTLIB:comdlg32.lib
        /NODEFAULTLIB:advapi32.lib"
)
