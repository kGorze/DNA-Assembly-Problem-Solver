#ifndef OPTYMALIZACJA_KOMBINATORYCZNA_LOGGER_TEST_H
#define OPTYMALIZACJA_KOMBINATORYCZNA_LOGGER_TEST_H

#include <gtest/gtest.h>
#include "utils/logging.h"

class LoggerTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Initialize logger for testing with DEBUG level
        Logger::init("test.log", LogLevel::DEBUG);
    }

    void TearDown() override {
        // Clean up logger after each test
        Logger::cleanup();
    }
};

#endif //OPTYMALIZACJA_KOMBINATORYCZNA_LOGGER_TEST_H 