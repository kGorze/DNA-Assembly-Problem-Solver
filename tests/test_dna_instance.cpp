#include <gtest/gtest.h>
#include "dna/dna_instance.h"

TEST(DNAInstanceTest, DefaultConstructor) {
    DNAInstance instance;
    EXPECT_EQ(instance.getN(), 0);
    EXPECT_EQ(instance.getK(), 0);
    EXPECT_EQ(instance.getDeltaK(), 0);
    EXPECT_EQ(instance.getLNeg(), 0);
    EXPECT_EQ(instance.getLPoz(), 0);
    EXPECT_FALSE(instance.isRepAllowed());
    EXPECT_EQ(instance.getProbablePositive(), 0.0);
    EXPECT_EQ(instance.getStartIndex(), 0);
    EXPECT_TRUE(instance.getSpectrum().empty());
}

TEST(DNAInstanceTest, ParameterizedConstructor) {
    DNAInstance instance(100, 10, 2, 5, 5, true, 0.8, 0);
    EXPECT_EQ(instance.getN(), 100);
    EXPECT_EQ(instance.getK(), 10);
    EXPECT_EQ(instance.getDeltaK(), 2);
    EXPECT_EQ(instance.getLNeg(), 5);
    EXPECT_EQ(instance.getLPoz(), 5);
    EXPECT_TRUE(instance.isRepAllowed());
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 0.8);
    EXPECT_EQ(instance.getStartIndex(), 0);
}

TEST(DNAInstanceTest, SettersAndGetters) {
    DNAInstance instance;
    
    instance.setN(50);
    EXPECT_EQ(instance.getN(), 50);
    
    instance.setK(8);
    EXPECT_EQ(instance.getK(), 8);
    
    instance.setDeltaK(1);
    EXPECT_EQ(instance.getDeltaK(), 1);
    
    instance.setLNeg(3);
    EXPECT_EQ(instance.getLNeg(), 3);
    
    instance.setLPoz(3);
    EXPECT_EQ(instance.getLPoz(), 3);
    
    instance.setRepAllowed(true);
    EXPECT_TRUE(instance.isRepAllowed());
    
    instance.setProbablePositive(0.7);
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 0.7);
    
    instance.setStartIndex(1);
    EXPECT_EQ(instance.getStartIndex(), 1);
}

TEST(DNAInstanceTest, SpectrumOperations) {
    DNAInstance instance;
    std::vector<std::string> spectrum = {"AAA", "AAT", "ATA", "TAA"};
    
    instance.setSpectrum(spectrum);
    EXPECT_EQ(instance.getSpectrum(), spectrum);
    EXPECT_EQ(instance.getSpectrum().size(), 4);
    
    instance.clearSpectrum();
    EXPECT_TRUE(instance.getSpectrum().empty());
}

TEST(DNAInstanceTest, InvalidParameters) {
    DNAInstance instance;
    
    // Test negative values
    instance.setN(-1);
    EXPECT_EQ(instance.getN(), 0);
    
    instance.setK(-1);
    EXPECT_EQ(instance.getK(), 0);
    
    instance.setDeltaK(-1);
    EXPECT_EQ(instance.getDeltaK(), 0);
    
    instance.setLNeg(-1);
    EXPECT_EQ(instance.getLNeg(), 0);
    
    instance.setLPoz(-1);
    EXPECT_EQ(instance.getLPoz(), 0);
    
    // Test probability bounds
    instance.setProbablePositive(-0.1);
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 0.0);
    
    instance.setProbablePositive(1.1);
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 1.0);
} 