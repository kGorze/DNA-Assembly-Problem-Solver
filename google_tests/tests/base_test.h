#pragma once

#include <gtest/gtest.h>
#include "utils/logging.h"
#include <iostream>

class BaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        const ::testing::TestInfo* const test_info = 
            ::testing::UnitTest::GetInstance()->current_test_info();
        std::cout << "\nStarting test: " << test_info->name() << "..." << std::endl;
        LOG_INFO("Starting test: {}", test_info->name());
    }

    void TearDown() override {
        const ::testing::TestInfo* const test_info = 
            ::testing::UnitTest::GetInstance()->current_test_info();
        std::cout << "Test " << test_info->name() << " completed." << std::endl;
        LOG_INFO("Test {} completed.", test_info->name());
    }
}; 