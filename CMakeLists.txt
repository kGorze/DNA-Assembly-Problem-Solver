cmake_minimum_required(VERSION 3.29)
project(optymalizacja_kombinatoryczna)

set(CMAKE_CXX_STANDARD 20)


add_executable(optymalizacja_kombinatoryczna src/main.cpp
        include/common.h
        src/generator/dna_generator.cpp
        include/generator/dna_generator.h
        src/naive/naive_reconstruction.cpp
        include/naive/naive_reconstruction.h
        src/benchmark/naive_benchmark.cpp
        include/benchmark/naive_benchmark.h
        src/metaheuristics/genetic_algorithm.cpp
        include/metaheuristics/genetic_algorithm.h
        src/metaheuristics/selection.cpp
        include/metaheuristics/selection.h
        include/metaheuristics/crossover.h
        src/metaheuristics/crossover.cpp
        include/metaheuristics/mutation.h
        src/metaheuristics/mutation.cpp
        include/metaheuristics/replacement.h
        src/metaheuristics/replacement.cpp
        include/metaheuristics/fitness.h
        src/metaheuristics/fitness.cpp
        src/metaheuristics/stopping_criteria.cpp
        include/metaheuristics/stopping_criteria.h
        include/metaheuristics/representation.h
        src/metaheuristics/representation.cpp
        include/benchmark/crossover_benchmark.h
        src/benchmark/crossover_benchmark.cpp
        src/utils/utility_functions.cpp
        include/utils/utility_functions.h
        include/configuration/genetic_algorithm_configuration.h
        src/configuration/genetic_algorithm_configuration.cpp
        include/utils/performance_profilling_framework.h
        src/utils/performance_profilling_framework.cpp
        include/metaheuristics/population_cache.h
        src/metaheuristics/population_cache.cpp
        include/utils/zobrist_hasher.h
)

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_SOURCE_DIR}/include)
