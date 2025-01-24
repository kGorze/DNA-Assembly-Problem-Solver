cmake_minimum_required(VERSION 3.29)
project(optymalizacja_kombinatoryczna)

# Set C++ standard and optimization flags
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set optimization level to O3 and enable static linking
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgcc -static-libstdc++ -static")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# Create library from source files
add_library(${PROJECT_NAME}_lib STATIC
        src/dna/dna_instance_builder.cpp
        src/dna/dna_instance_io.cpp
        src/generator/dna_generator.cpp
        src/naive/naive_reconstruction.cpp
        src/utils/utility_functions.cpp
        src/utils/performance_profilling_framework.cpp
        src/utils/logging.cpp
        src/metaheuristics/genetic_algorithm.cpp
        src/metaheuristics/genetic_algorithm_runner.cpp
        src/metaheuristics/crossover.cpp
        src/metaheuristics/mutation.cpp
        src/metaheuristics/selection.cpp
        src/metaheuristics/replacement.cpp
        src/metaheuristics/fitness.cpp
        src/metaheuristics/stopping_criteria.cpp
        src/metaheuristics/population_cache.cpp
        src/metaheuristics/adaptive_crossover.cpp
        src/metaheuristics/path_analyzer.cpp
        src/metaheuristics/representation.cpp
        src/configuration/genetic_algorithm_configuration.cpp
        src/benchmark/crossover_benchmark.cpp
        src/benchmark/naive_benchmark.cpp
        src/benchmark/adaptive_crossover_benchmark.cpp
        src/tuning/parameter_tuning_manager.cpp
        src/tuning/parameter_parser.cpp
        src/tuning/tuning_structures.cpp
)

target_include_directories(${PROJECT_NAME}_lib PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
)

# Main executable
add_executable(${PROJECT_NAME}
        src/main.cpp
)

target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_lib)

# Google Test setup
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(googletest)

enable_testing()
add_subdirectory(google_tests)