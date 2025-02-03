// tests/test_main.cpp
#include <gtest/gtest.h>
#include "utils/logging.h"
#include <filesystem>
#include "base_test.h"

class TestEnvironment : public ::testing::Environment {
public:
    ~TestEnvironment() override = default;

    void SetUp() override {
        // Clean up any existing test log file
        if (std::filesystem::exists("test.log")) {
            std::filesystem::remove("test.log");
        }
        // Initialize logger for tests
        Logger::initialize("test.log");
        Logger::setLogLevel(LogLevel::DEBUG);
    }

    void TearDown() override {
        // Clean up logger after tests
        Logger::cleanup();
        // Remove the test log file
        if (std::filesystem::exists("test.log")) {
            std::filesystem::remove("test.log");
        }
    }
};

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    testing::AddGlobalTestEnvironment(new TestEnvironment);
    return RUN_ALL_TESTS();
}