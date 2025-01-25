#include <gtest/gtest.h>
#include "generator/dna_generator.h"
#include "dna/dna_instance.h"

TEST(DNAGeneratorTest, GenerateDNA) {
    DNAGenerator generator;
    
    // Test with valid length
    std::string dna = generator.generateDNA(10);
    EXPECT_EQ(dna.length(), 10);
    
    // Verify that only valid nucleotides are used
    for (char c : dna) {
        EXPECT_TRUE(c == 'A' || c == 'T' || c == 'G' || c == 'C');
    }
    
    // Test with zero length
    EXPECT_TRUE(generator.generateDNA(0).empty());
    
    // Test with negative length
    EXPECT_TRUE(generator.generateDNA(-1).empty());
}

TEST(DNAGeneratorTest, GenerateRandomInstance) {
    DNAGenerator generator;
    
    // Test with valid parameters
    auto instance = generator.generateRandomInstance(100, 10, 2, 5, 5, true, 0.8);
    EXPECT_EQ(instance.getN(), 100);
    EXPECT_EQ(instance.getK(), 10);
    EXPECT_EQ(instance.getDeltaK(), 2);
    EXPECT_EQ(instance.getLNeg(), 5);
    EXPECT_EQ(instance.getLPoz(), 5);
    EXPECT_TRUE(instance.isRepAllowed());
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 0.8);
    EXPECT_FALSE(instance.getSpectrum().empty());
    
    // Test with invalid parameters
    auto invalidInstance = generator.generateRandomInstance(-1, -1, -1, -1, -1, false, -0.1);
    EXPECT_EQ(invalidInstance.getN(), 0);
    EXPECT_EQ(invalidInstance.getK(), 0);
    EXPECT_EQ(invalidInstance.getDeltaK(), 0);
    EXPECT_EQ(invalidInstance.getLNeg(), 0);
    EXPECT_EQ(invalidInstance.getLPoz(), 0);
    EXPECT_FALSE(invalidInstance.isRepAllowed());
    EXPECT_DOUBLE_EQ(invalidInstance.getProbablePositive(), 0.0);
    EXPECT_TRUE(invalidInstance.getSpectrum().empty());
}

TEST(DNAGeneratorTest, SaveToFile) {
    DNAGenerator generator;
    auto instance = generator.generateRandomInstance(100, 10, 2, 5, 5, true, 0.8);
    
    // Test with valid filename
    EXPECT_TRUE(generator.saveToFile(instance, "test_instance.txt"));
    
    // Test with empty filename
    EXPECT_FALSE(generator.saveToFile(instance, ""));
    
    // Clean up
    std::remove("test_instance.txt");
} 