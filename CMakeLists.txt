cmake_minimum_required(VERSION 3.22)
project(optymalizacja_kombinatoryczna)

# Set C++ standard and optimization flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set optimization level to O3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Find Python (if needed)
find_package(Python3 COMPONENTS Interpreter REQUIRED)

# Define include directories at the top level
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}
)

# Source files
set(SOURCES
    src/configuration/genetic_algorithm_configuration.cpp
    src/metaheuristics/mutation_impl.cpp
    src/metaheuristics/mutation.cpp
    src/metaheuristics/selection.cpp
    src/metaheuristics/crossover_impl.cpp
    src/metaheuristics/crossover.cpp
    src/metaheuristics/permutation_representation.cpp
    src/metaheuristics/genetic_algorithm.cpp
    src/metaheuristics/genetic_algorithm_runner.cpp
    src/metaheuristics/adaptive_crossover.cpp
    src/metaheuristics/fitness_impl.cpp
    src/metaheuristics/fitness.cpp
    src/metaheuristics/individual.cpp
    src/metaheuristics/stopping_criteria.cpp
    src/metaheuristics/path_analyzer.cpp
    src/dna/dna_instance.cpp
    src/dna/dna_instance_io.cpp
    src/dna/dna_instance_builder.cpp
    src/dna/error_introduction.cpp
    src/dna/dna_generator.cpp
    src/dna/spectrum_generator.cpp
    src/utils/logging.cpp
    src/utils/utility_functions.cpp
    src/utils/performance_profilling_framework.cpp
    src/tuning/parameter_tuning_manager.cpp
    src/tuning/parameter_parser.cpp
    src/tuning/tuning_structures.cpp
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

# Copy configuration files to build directory
file(COPY ${CMAKE_SOURCE_DIR}/config DESTINATION ${CMAKE_BINARY_DIR})

# Set working directory for the executable
set_target_properties(optymalizacja_kombinatoryczna PROPERTIES
    VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
)

# Only include tests if ENABLE_TESTING is set
option(ENABLE_TESTING "Enable testing" OFF)
if(ENABLE_TESTING)
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
endif()