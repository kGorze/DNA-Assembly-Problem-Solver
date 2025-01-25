#include <gtest/gtest.h>
#include <mutex>
#include <memory>
#include "generator/dna_generator.h"
#include "utils/logging.h"

namespace {
    std::mutex test_mutex;
}

class DNAGeneratorTest : public ::testing::Test {
protected:
    void SetUp() override {
        generator = std::make_unique<DNAGenerator>();
        Logger::initialize("test.log");
    }

    void TearDown() override {
        generator.reset();
        Logger::cleanup();
    }

    std::unique_ptr<DNAGenerator> generator;
};

TEST_F(DNAGeneratorTest, GeneratesCorrectLength) {
    std::lock_guard<std::mutex> lock(test_mutex);
    const int n = 100;
    try {
        std::string dna = generator->generateDNA(n, true);
        EXPECT_EQ(dna.length(), n) << "Generated DNA length does not match expected length";
    } catch (const std::exception& e) {
        FAIL() << "Exception thrown during DNA generation: " << e.what();
    }
}

TEST_F(DNAGeneratorTest, ContainsValidNucleotides) {
    std::lock_guard<std::mutex> lock(test_mutex);
    try {
        std::string dna = generator->generateDNA(50, true);
        for (char c : dna) {
            EXPECT_TRUE(c == 'A' || c == 'C' || c == 'G' || c == 'T')
                << "Invalid nucleotide found: " << c;
        }
    } catch (const std::exception& e) {
        FAIL() << "Exception thrown during DNA generation: " << e.what();
    }
}

TEST(DNAInstanceTest, BasicSettersAndGetters) {
    std::lock_guard<std::mutex> lock(test_mutex);
    try {
        DNAInstance instance;
        instance.setN(100);
        instance.setK(10);
        instance.setDeltaK(2);
        
        EXPECT_EQ(instance.getN(), 100) << "N value mismatch";
        EXPECT_EQ(instance.getK(), 10) << "K value mismatch";
        EXPECT_EQ(instance.getDeltaK(), 2) << "DeltaK value mismatch";
    } catch (const std::exception& e) {
        FAIL() << "Exception thrown during DNAInstance test: " << e.what();
    }
}

TEST(SpectrumGeneratorTest, GeneratesValidSpectrum) {
    std::lock_guard<std::mutex> lock(test_mutex);
    try {
        SpectrumGenerator generator;
        const std::string dna = "ACGTACGT";
        const int k = 3;
        const int deltaK = 1;
        
        auto spectrum = generator.generateSpectrum(dna, k, deltaK);
        EXPECT_FALSE(spectrum.empty()) << "Generated spectrum is empty";
        
        for (const auto& oligo : spectrum) {
            EXPECT_GE(oligo.length(), k - deltaK) 
                << "Oligo length " << oligo.length() << " is less than minimum length " << (k - deltaK);
            EXPECT_LE(oligo.length(), k + deltaK)
                << "Oligo length " << oligo.length() << " is greater than maximum length " << (k + deltaK);
        }
    } catch (const std::exception& e) {
        FAIL() << "Exception thrown during spectrum generation: " << e.what();
    }
}