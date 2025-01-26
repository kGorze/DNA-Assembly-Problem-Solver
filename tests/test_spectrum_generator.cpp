#include <gtest/gtest.h>
#include "../include/generator/dna_generator.h"
#include "../include/dna/dna_instance.h"

class SpectrumGeneratorTest : public ::testing::Test {
protected:
    void SetUp() override {}
    void TearDown() override {}
};

TEST_F(SpectrumGeneratorTest, GenerateSpectrum_Test) {
    DNAInstance instance;
    instance.setK(3);
    instance.setDNA("ACGTACGT");
    
    DNAGenerator generator;
    auto spectrum = generator.generateDNASpectrum(instance);
    
    ASSERT_FALSE(spectrum.empty());
    for (const auto& kmer : spectrum) {
        ASSERT_EQ(kmer.length(), instance.getK());
    }
}

TEST_F(SpectrumGeneratorTest, GenerateSpectrumWithRepetitions_Test) {
    DNAInstance instance;
    instance.setK(3);
    instance.setDNA("AAAAAAA");
    
    DNAGenerator generator;
    auto spectrum = generator.generateDNASpectrum(instance);
    
    ASSERT_FALSE(spectrum.empty());
    for (const auto& kmer : spectrum) {
        ASSERT_EQ(kmer.length(), instance.getK());
    }
}

TEST_F(SpectrumGeneratorTest, GenerateSpectrumWithoutRepetitions_Test) {
    DNAInstance instance;
    instance.setK(3);
    instance.setDNA("ACGTACGT");
    
    DNAGenerator generator;
    auto spectrum = generator.generateDNASpectrum(instance);
    
    ASSERT_FALSE(spectrum.empty());
    for (const auto& kmer : spectrum) {
        ASSERT_EQ(kmer.length(), instance.getK());
    }
}

TEST_F(SpectrumGeneratorTest, GenerateSpectrumWithDifferentLengths_Test) {
    DNAInstance instance;
    instance.setK(4);
    instance.setDNA("ACGTACGT");
    
    DNAGenerator generator;
    auto spectrum = generator.generateDNASpectrum(instance);
    
    ASSERT_FALSE(spectrum.empty());
    for (const auto& kmer : spectrum) {
        ASSERT_EQ(kmer.length(), instance.getK());
    }
} 