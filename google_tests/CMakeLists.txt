add_executable(Google_Tests_run
        tests/test_main.cpp
        tests/your_tests.cpp
)

target_link_libraries(Google_Tests_run PRIVATE
        gtest
        gtest_main
        ${PROJECT_NAME}_lib
)

target_include_directories(Google_Tests_run PRIVATE ${CMAKE_SOURCE_DIR}/include)

add_test(NAME GoogleTests COMMAND Google_Tests_run)