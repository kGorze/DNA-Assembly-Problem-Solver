#include <gtest/gtest.h>
#include "generator/spectrum_generator.h"
#include "dna/dna_instance.h"

TEST(SpectrumGeneratorTest, GenerateSpectrum) {
    SpectrumGenerator generator;
    DNAInstance instance(100, 10, 2, 5, 5, true, 0.8, 0);
    
    // Test with valid parameters
    auto spectrum = generator.generateSpectrum(instance);
    EXPECT_FALSE(spectrum.empty());
    
    // Verify that all oligonucleotides have valid length
    for (const auto& oligo : spectrum) {
        EXPECT_GE(oligo.length(), instance.getK() - instance.getDeltaK());
        EXPECT_LE(oligo.length(), instance.getK() + instance.getDeltaK());
        
        // Verify that only valid nucleotides are used
        for (char c : oligo) {
            EXPECT_TRUE(c == 'A' || c == 'T' || c == 'G' || c == 'C');
        }
    }
    
    // Test with invalid parameters
    DNAInstance invalidInstance;
    auto emptySpectrum = generator.generateSpectrum(invalidInstance);
    EXPECT_TRUE(emptySpectrum.empty());
}

TEST(SpectrumGeneratorTest, GenerateSpectrumWithRepetitions) {
    SpectrumGenerator generator;
    DNAInstance instance(100, 10, 2, 5, 5, true, 0.8, 0);
    instance.setRepAllowed(true);
    
    auto spectrum = generator.generateSpectrum(instance);
    EXPECT_FALSE(spectrum.empty());
    
    // Count occurrences of each oligonucleotide
    std::map<std::string, int> counts;
    for (const auto& oligo : spectrum) {
        counts[oligo]++;
    }
    
    // Check if there are any repetitions
    bool hasRepetitions = false;
    for (const auto& [oligo, count] : counts) {
        if (count > 1) {
            hasRepetitions = true;
            break;
        }
    }
    EXPECT_TRUE(hasRepetitions);
}

TEST(SpectrumGeneratorTest, GenerateSpectrumWithoutRepetitions) {
    SpectrumGenerator generator;
    DNAInstance instance(100, 10, 2, 5, 5, false, 0.8, 0);
    instance.setRepAllowed(false);
    
    auto spectrum = generator.generateSpectrum(instance);
    EXPECT_FALSE(spectrum.empty());
    
    // Count occurrences of each oligonucleotide
    std::map<std::string, int> counts;
    for (const auto& oligo : spectrum) {
        counts[oligo]++;
    }
    
    // Verify that there are no repetitions
    for (const auto& [oligo, count] : counts) {
        EXPECT_EQ(count, 1);
    }
}

TEST(SpectrumGeneratorTest, GenerateSpectrumWithDifferentLengths) {
    SpectrumGenerator generator;
    DNAInstance instance(100, 10, 2, 5, 5, false, 0.8, 0);
    
    auto spectrum = generator.generateSpectrum(instance);
    EXPECT_FALSE(spectrum.empty());
    
    // Count occurrences of each length
    std::map<size_t, int> lengthCounts;
    for (const auto& oligo : spectrum) {
        lengthCounts[oligo.length()]++;
    }
    
    // Verify that we have oligonucleotides of different lengths
    EXPECT_GT(lengthCounts.size(), 1);
    
    // Verify that all lengths are within the allowed range
    for (const auto& [length, count] : lengthCounts) {
        EXPECT_GE(length, instance.getK() - instance.getDeltaK());
        EXPECT_LE(length, instance.getK() + instance.getDeltaK());
    }
} 