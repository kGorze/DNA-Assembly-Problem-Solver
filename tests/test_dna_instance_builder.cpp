#include <gtest/gtest.h>
#include "dna/dna_instance_builder.h"
#include "generator/dna_generator.h"
#include "generator/spectrum_generator.h"
#include "dna/error_introduction.h"

TEST(DNAInstanceBuilderTest, BuildWithValidParameters) {
    DNAInstanceBuilder builder;
    
    // Set valid parameters
    builder.setN(100)
           .setK(10)
           .setDeltaK(2)
           .setLNeg(5)
           .setLPoz(5)
           .setRepAllowed(true)
           .setProbablePositive(0.8)
           .setStartIndex(0);
    
    // Build instance
    auto instance = builder.build();
    
    // Verify parameters
    EXPECT_EQ(instance.getN(), 100);
    EXPECT_EQ(instance.getK(), 10);
    EXPECT_EQ(instance.getDeltaK(), 2);
    EXPECT_EQ(instance.getLNeg(), 5);
    EXPECT_EQ(instance.getLPoz(), 5);
    EXPECT_TRUE(instance.isRepAllowed());
    EXPECT_EQ(instance.getProbablePositive(), 0.8);
    EXPECT_EQ(instance.getStartIndex(), 0);
}

TEST(DNAInstanceBuilderTest, BuildWithInvalidParameters) {
    DNAInstanceBuilder builder;
    
    // Set invalid parameters
    builder.setN(-1)
           .setK(-1)
           .setDeltaK(-1)
           .setLNeg(-1)
           .setLPoz(-1)
           .setRepAllowed(false)
           .setProbablePositive(-0.1)
           .setStartIndex(-1);
    
    // Build instance should throw
    EXPECT_THROW(builder.build(), std::invalid_argument);
}

TEST(DNAInstanceBuilderTest, BuildWithDefaultParameters) {
    DNAInstanceBuilder builder;
    
    // Build instance with default parameters
    auto instance = builder.build();
    
    // Verify default parameters
    EXPECT_EQ(instance.getN(), 0);
    EXPECT_EQ(instance.getK(), 0);
    EXPECT_EQ(instance.getDeltaK(), 0);
    EXPECT_EQ(instance.getLNeg(), 0);
    EXPECT_EQ(instance.getLPoz(), 0);
    EXPECT_FALSE(instance.isRepAllowed());
    EXPECT_EQ(instance.getProbablePositive(), 0.0);
    EXPECT_EQ(instance.getStartIndex(), 0);
}

TEST(DNAInstanceBuilderTest, BuildWithDNA) {
    DNAInstanceBuilder builder;
    
    // Set valid parameters
    builder.setN(12)
           .setK(4)
           .setDeltaK(1)
           .setLNeg(2)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.8)
           .setStartIndex(0);
    
    // Set DNA
    builder.setDNA("ATCGATCGATCG");
    
    // Build instance
    auto instance = builder.build();
    
    // Verify DNA
    EXPECT_EQ(instance.getDNA(), "ATCGATCGATCG");
}

TEST(DNAInstanceBuilderTest, BuildWithSpectrum) {
    DNAInstanceBuilder builder;
    
    // Set valid parameters
    builder.setN(12)
           .setK(4)
           .setDeltaK(1)
           .setLNeg(2)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.8)
           .setStartIndex(0);
    
    // Set spectrum
    std::vector<std::string> spectrum = {"ATCG", "TCGA", "CGAT"};
    builder.setSpectrum(spectrum);
    
    // Build instance
    auto instance = builder.build();
    
    // Verify spectrum
    EXPECT_EQ(instance.getSpectrum(), spectrum);
}

TEST(DNAInstanceBuilderTest, BuildWithInvalidDNA) {
    DNAInstanceBuilder builder;
    
    // Set valid parameters
    builder.setN(12)
           .setK(4)
           .setDeltaK(1)
           .setLNeg(2)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.8)
           .setStartIndex(0);
    
    // Set invalid DNA
    builder.setDNA("INVALID123");
    
    // Build instance should throw
    EXPECT_THROW(builder.build(), std::invalid_argument);
}

TEST(DNAInstanceBuilderTest, BuildWithInvalidSpectrum) {
    DNAInstanceBuilder builder;
    
    // Set valid parameters
    builder.setN(12)
           .setK(4)
           .setDeltaK(1)
           .setLNeg(2)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.8)
           .setStartIndex(0);
    
    // Set invalid spectrum
    std::vector<std::string> spectrum = {"INVALID", "123", "XYZ"};
    builder.setSpectrum(spectrum);
    
    // Build instance should throw
    EXPECT_THROW(builder.build(), std::invalid_argument);
} 