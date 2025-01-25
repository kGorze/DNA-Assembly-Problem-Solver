set(TEST_SOURCES
    test_dna_instance.cpp
    test_dna_generator.cpp
    test_spectrum_generator.cpp
    test_error_introduction.cpp
    test_genetic_algorithm.cpp
)

add_executable(optymalizacja_kombinatoryczna_test ${TEST_SOURCES})

target_link_libraries(optymalizacja_kombinatoryczna_test PRIVATE
    optymalizacja_kombinatoryczna_lib
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(optymalizacja_kombinatoryczna_test) 