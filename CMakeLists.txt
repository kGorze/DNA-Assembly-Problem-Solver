cmake_minimum_required(VERSION 3.16)
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

# Create library from source files
add_library(${PROJECT_NAME}_lib STATIC
        src/dna/dna_instance.cpp
        src/dna/dna_instance_builder.cpp
        src/dna/error_introduction.cpp
        src/dna/dna_instance_io.cpp
        src/generator/dna_generator.cpp
        src/generator/spectrum_generator.cpp
        src/utils/logging.cpp
        src/metaheuristics/genetic_algorithm.cpp
        src/metaheuristics/representation.cpp
        src/metaheuristics/population_cache.cpp
)

# Set include directories for the library
target_include_directories(${PROJECT_NAME}_lib PUBLIC 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src
)

# Main executable
add_executable(${PROJECT_NAME}
        src/main.cpp
)

# Link the library
target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}_lib)

# Google Test setup
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

# Add the Google Tests executable
add_executable(${PROJECT_NAME}_test
    google_tests/tests/test_main.cpp
    google_tests/tests/your_tests.cpp
)

target_link_libraries(${PROJECT_NAME}_test PRIVATE
    ${PROJECT_NAME}_lib
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME}_test)