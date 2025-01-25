cmake_minimum_required(VERSION 3.22)
project(optymalizacja_kombinatoryczna)

# Set C++ standard and optimization flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set optimization level to O3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Define include directories at the top level
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}
)

# Source files
set(SOURCES
    src/dna/dna_instance.cpp
    src/dna/dna_instance_builder.cpp
    src/dna/error_introduction.cpp
    src/dna/dna_instance_io.cpp
    src/generator/dna_generator.cpp
    src/generator/spectrum_generator.cpp
    src/utils/logging.cpp
    src/metaheuristics/genetic_algorithm.cpp
    src/metaheuristics/population_cache.cpp
)

# Create library
add_library(optymalizacja_kombinatoryczna_lib STATIC ${SOURCES})

# Set include directories for the library
target_include_directories(optymalizacja_kombinatoryczna_lib PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

# Main executable
add_executable(optymalizacja_kombinatoryczna
        src/main.cpp
)

# Link the library
target_link_libraries(optymalizacja_kombinatoryczna PRIVATE optymalizacja_kombinatoryczna_lib)

# Google Test setup
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    DOWNLOAD_EXTRACT_TIMESTAMP true
)
FetchContent_MakeAvailable(googletest)

# Enable testing
enable_testing()

# Add tests
add_subdirectory(tests)