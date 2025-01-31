// tests/test_main.cpp
#include <gtest/gtest.h>
#include "utils/logging.h"

class TestEnvironment : public ::testing::Environment {
public:
    ~TestEnvironment() override = default;

    void SetUp() override {
        // Inicjalizacja loggera dla testów
        Logger::initialize("test.log");
        Logger::setLogLevel(LogLevel::DEBUG);
    }

    void TearDown() override {
        // Czyszczenie loggera po testach
        Logger::cleanup();
    }
};

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    testing::AddGlobalTestEnvironment(new TestEnvironment);
    return RUN_ALL_TESTS();
}