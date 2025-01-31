#include <gtest/gtest.h>
#include "dna/dna_instance.h"

TEST(DNAInstanceTest, CreateEmptyInstance) {
    DNAInstance instance;
    EXPECT_EQ(instance.getSpectrum().size(), 0);
}

TEST(DNAInstanceTest, AddDNASequence) {
    DNAInstance instance;
    instance.setDNA("ACGT");
    EXPECT_EQ(instance.getDNA(), "ACGT");
}

TEST(DNAInstanceTest, SetAndGetSpectrum) {
    DNAInstance instance;
    std::vector<std::string> spectrum = {"ACGT", "CGTA"};
    instance.setSpectrum(spectrum);
    EXPECT_EQ(instance.getSpectrum(), spectrum);
}