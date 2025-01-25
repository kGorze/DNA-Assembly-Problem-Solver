#include <gtest/gtest.h>
#include <memory>
#include "utils/logging.h"

class GlobalEnvironment : public testing::Environment {
public:
    void SetUp() override {
        Logger::initialize("test.log");
    }

    void TearDown() override {
        Logger::cleanup();
    }
};

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    testing::AddGlobalTestEnvironment(new GlobalEnvironment);
    const int result = RUN_ALL_TESTS();
    return result;
}