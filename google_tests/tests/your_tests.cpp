#include <gtest/gtest.h>
#include "generator/dna_generator.h"

class DNAGeneratorTest : public ::testing::Test {
protected:
    DNAGenerator generator;
};

TEST_F(DNAGeneratorTest, GeneratesCorrectLength) {
    int n = 100;
    std::string dna = generator.generateDNA(n, true);
    EXPECT_EQ(dna.length(), n);
}

TEST_F(DNAGeneratorTest, ContainsValidNucleotides) {
    std::string dna = generator.generateDNA(50, true);
    for (char c : dna) {
        EXPECT_TRUE(c == 'A' || c == 'C' || c == 'G' || c == 'T');
    }
}

TEST(DNAInstanceTest, BasicSettersAndGetters) {
    DNAInstance instance;
    instance.setN(100);
    instance.setK(10);
    instance.setDeltaK(2);
    
    EXPECT_EQ(instance.getN(), 100);
    EXPECT_EQ(instance.getK(), 10);
    EXPECT_EQ(instance.getDeltaK(), 2);
}

TEST(SpectrumGeneratorTest, GeneratesValidSpectrum) {
    SpectrumGenerator generator;
    std::string dna = "ACGTACGT";
    int k = 3;
    int deltaK = 1;
    
    auto spectrum = generator.generateSpectrum(dna, k, deltaK);
    EXPECT_FALSE(spectrum.empty());
    
    for (const auto& oligo : spectrum) {
        EXPECT_GE(oligo.length(), k - deltaK);
        EXPECT_LE(oligo.length(), k + deltaK);
    }
}