#include <gtest/gtest.h>
#include "dna/error_introduction.h"
#include "dna/dna_instance.h"
#include <algorithm>

TEST(ErrorIntroductionTest, PositiveErrorIntroducer) {
    PositiveErrorIntroducer introducer;
    DNAInstance instance(100, 10, 2, 5, 5, true, 0.8, 0);
    
    // Create a spectrum
    std::vector<std::string> spectrum = {
        "AAAA", "AAAT", "AATA", "ATAA", "TAAA"
    };
    instance.setSpectrum(spectrum);
    
    // Introduce positive errors
    introducer.introduceErrors(instance);
    const auto& modifiedSpectrum = instance.getSpectrum();
    
    // Verify that the modified spectrum is not empty
    EXPECT_FALSE(modifiedSpectrum.empty());
    
    // Verify that the modified spectrum has more elements than the original
    EXPECT_GE(modifiedSpectrum.size(), spectrum.size());
    
    // Verify that all original elements are still present
    for (const auto& oligo : spectrum) {
        EXPECT_TRUE(std::find(modifiedSpectrum.begin(), modifiedSpectrum.end(), oligo) != modifiedSpectrum.end());
    }
    
    // Verify that added oligonucleotides have valid length
    for (const auto& oligo : modifiedSpectrum) {
        EXPECT_GE(oligo.length(), instance.getK() - instance.getDeltaK());
        EXPECT_LE(oligo.length(), instance.getK() + instance.getDeltaK());
        
        // Verify that only valid nucleotides are used
        for (char c : oligo) {
            EXPECT_TRUE(c == 'A' || c == 'T' || c == 'G' || c == 'C');
        }
    }
}

TEST(ErrorIntroductionTest, NegativeErrorIntroducer) {
    NegativeErrorIntroducer introducer;
    DNAInstance instance(100, 10, 2, 5, 5, true, 0.8, 0);
    
    // Create a spectrum
    std::vector<std::string> spectrum = {
        "AAAA", "AAAT", "AATA", "ATAA", "TAAA",
        "TTTT", "TTTA", "TTAT", "TATT", "ATTT"
    };
    instance.setSpectrum(spectrum);
    
    // Introduce negative errors
    introducer.introduceErrors(instance);
    const auto& modifiedSpectrum = instance.getSpectrum();
    
    // Verify that the modified spectrum is not empty
    EXPECT_FALSE(modifiedSpectrum.empty());
    
    // Verify that the modified spectrum has fewer elements than the original
    EXPECT_LE(modifiedSpectrum.size(), spectrum.size());
    
    // Verify that removed oligonucleotides are actually missing
    bool hasRemovals = false;
    for (const auto& oligo : spectrum) {
        if (std::find(modifiedSpectrum.begin(), modifiedSpectrum.end(), oligo) == modifiedSpectrum.end()) {
            hasRemovals = true;
            break;
        }
    }
    EXPECT_TRUE(hasRemovals);
    
    // Verify that remaining oligonucleotides are from the original spectrum
    for (const auto& oligo : modifiedSpectrum) {
        EXPECT_TRUE(std::find(spectrum.begin(), spectrum.end(), oligo) != spectrum.end());
    }
}

TEST(ErrorIntroductionTest, EmptySpectrum) {
    PositiveErrorIntroducer positiveIntroducer;
    NegativeErrorIntroducer negativeIntroducer;
    DNAInstance instance;
    
    // Test with empty spectrum
    positiveIntroducer.introduceErrors(instance);
    const auto& modifiedSpectrumPos = instance.getSpectrum();
    EXPECT_TRUE(modifiedSpectrumPos.empty());
    
    negativeIntroducer.introduceErrors(instance);
    const auto& modifiedSpectrumNeg = instance.getSpectrum();
    EXPECT_TRUE(modifiedSpectrumNeg.empty());
}

TEST(ErrorIntroductionTest, InvalidParameters) {
    PositiveErrorIntroducer positiveIntroducer;
    NegativeErrorIntroducer negativeIntroducer;
    
    // Create instance with invalid parameters
    DNAInstance instance(-1, -1, -1, -1, -1, false, -0.1, -1);
    std::vector<std::string> spectrum = {"AAAA", "TTTT"};
    instance.setSpectrum(spectrum);
    
    // Test error introduction with invalid parameters
    positiveIntroducer.introduceErrors(instance);
    const auto& modifiedSpectrumPos = instance.getSpectrum();
    EXPECT_EQ(modifiedSpectrumPos, spectrum);  // Should return original spectrum
    
    negativeIntroducer.introduceErrors(instance);
    const auto& modifiedSpectrumNeg = instance.getSpectrum();
    EXPECT_EQ(modifiedSpectrumNeg, spectrum);  // Should return original spectrum
} 