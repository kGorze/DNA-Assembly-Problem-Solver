#include <gtest/gtest.h>
#include "dna/dna_instance_io.h"
#include "dna/dna_instance.h"
#include <filesystem>
#include <fstream>

class DNAInstanceIOTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a test instance
        instance = std::make_unique<DNAInstance>(100, 10, 2, 5, 5, true, 0.8, 0);
        instance->setDNA("ATCGATCGATCG");
        std::vector<std::string> spectrum = {"ATCG", "TCGA", "CGAT"};
        instance->setSpectrum(spectrum);
        
        // Set up test file path
        testFilePath = "test_instance.txt";
    }
    
    void TearDown() override {
        // Clean up test file if it exists
        if (std::filesystem::exists(testFilePath)) {
            std::filesystem::remove(testFilePath);
        }
    }
    
    std::unique_ptr<DNAInstance> instance;
    std::string testFilePath;
};

TEST_F(DNAInstanceIOTest, SaveAndLoadInstance) {
    // Save the instance
    EXPECT_NO_THROW(InstanceIO::saveInstance(*instance, testFilePath));
    
    // Verify file exists
    EXPECT_TRUE(std::filesystem::exists(testFilePath));
    
    // Load the instance
    auto loadedInstance = std::make_unique<DNAInstance>();
    EXPECT_NO_THROW(InstanceIO::loadInstance(*loadedInstance, testFilePath));
    
    // Verify loaded instance matches original
    EXPECT_EQ(loadedInstance->getN(), instance->getN());
    EXPECT_EQ(loadedInstance->getK(), instance->getK());
    EXPECT_EQ(loadedInstance->getDeltaK(), instance->getDeltaK());
    EXPECT_EQ(loadedInstance->getLNeg(), instance->getLNeg());
    EXPECT_EQ(loadedInstance->getLPoz(), instance->getLPoz());
    EXPECT_EQ(loadedInstance->isRepAllowed(), instance->isRepAllowed());
    EXPECT_EQ(loadedInstance->getProbablePositive(), instance->getProbablePositive());
    EXPECT_EQ(loadedInstance->getStartIndex(), instance->getStartIndex());
    EXPECT_EQ(loadedInstance->getDNA(), instance->getDNA());
    EXPECT_EQ(loadedInstance->getSpectrum(), instance->getSpectrum());
}

TEST_F(DNAInstanceIOTest, SaveToNonexistentDirectory) {
    std::string nonexistentPath = "nonexistent/directory/instance.txt";
    EXPECT_THROW(InstanceIO::saveInstance(*instance, nonexistentPath), std::runtime_error);
}

TEST_F(DNAInstanceIOTest, LoadNonexistentFile) {
    std::string nonexistentFile = "nonexistent_file.txt";
    auto loadedInstance = std::make_unique<DNAInstance>();
    EXPECT_THROW(InstanceIO::loadInstance(*loadedInstance, nonexistentFile), std::runtime_error);
}

TEST_F(DNAInstanceIOTest, SaveEmptyInstance) {
    auto emptyInstance = std::make_unique<DNAInstance>();
    EXPECT_NO_THROW(InstanceIO::saveInstance(*emptyInstance, testFilePath));
    
    // Load and verify empty instance
    auto loadedInstance = std::make_unique<DNAInstance>();
    EXPECT_NO_THROW(InstanceIO::loadInstance(*loadedInstance, testFilePath));
    
    EXPECT_TRUE(loadedInstance->getDNA().empty());
    EXPECT_TRUE(loadedInstance->getSpectrum().empty());
}

TEST_F(DNAInstanceIOTest, LoadCorruptedFile) {
    // Create a corrupted file
    std::ofstream corruptedFile(testFilePath);
    corruptedFile << "This is not a valid DNA instance file format";
    corruptedFile.close();
    
    auto loadedInstance = std::make_unique<DNAInstance>();
    EXPECT_THROW(InstanceIO::loadInstance(*loadedInstance, testFilePath), std::runtime_error);
} 