add_executable(optymalizacja_kombinatoryczna_test
        tests/test_main.cpp
        tests/dna_tests/dna_instance_test.cpp
        tests/metaheuristics/adaptive_crossover_test.cpp
        tests/dna_tests/dna_instance_builder_test.cpp
        tests/dna_tests/dna_generator_test.cpp
        tests/utils_tests/random_test.cpp
        tests/dna_tests/dna_instance_io_test.cpp
        tests/dna_tests/error_introduction_test.cpp
        tests/dna_tests/spectrum_generator_test.cpp
)

target_link_libraries(optymalizacja_kombinatoryczna_test PRIVATE
        GTest::gtest
        GTest::gtest_main
        GTest::gmock
        GTest::gmock_main
        ${PROJECT_NAME}_lib
)

target_include_directories(optymalizacja_kombinatoryczna_test PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src
)

include(GoogleTest)
gtest_discover_tests(optymalizacja_kombinatoryczna_test)