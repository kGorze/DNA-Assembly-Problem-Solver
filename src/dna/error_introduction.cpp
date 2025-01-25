#include "dna/error_introduction.h"
#include "utils/logging.h"
#include <algorithm>
#include <random>

void NegativeErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (!validateInstance(instance)) {
        return;
    }

    try {
        auto spectrum = instance.getSpectrum();
        if (spectrum.size() <= m_lNeg) {
            LOG_WARNING("Cannot remove " + std::to_string(m_lNeg) +
                       " k-mers from spectrum of size " + std::to_string(spectrum.size()));
            return;
        }

        // Randomly remove k-mers
        std::shuffle(spectrum.begin(), spectrum.end(), m_rng);
        spectrum.erase(spectrum.end() - m_lNeg, spectrum.end());
        instance.setSpectrum(spectrum);
    } catch (const std::exception& e) {
        LOG_ERROR("Error in negative error introduction: " + std::string(e.what()));
    }
}

void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (!validateInstance(instance)) {
        return;
    }

    try {
        auto spectrum = instance.getSpectrum();
        const int deltaK = instance.getDeltaK();
        
        // Add random k-mers
        for (int i = 0; i < m_lPoz; ++i) {
            int len = m_k + (m_rng() % (2 * deltaK + 1)) - deltaK;
            std::string newKmer = generateRandomKmer(len);
            spectrum.push_back(newKmer);
        }

        std::sort(spectrum.begin(), spectrum.end());
        instance.setSpectrum(spectrum);
    } catch (const std::exception& e) {
        LOG_ERROR("Error in positive error introduction: " + std::string(e.what()));
    }
}

std::string PositiveErrorIntroducer::generateRandomKmer(int length) const {
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::string kmer;
    kmer.reserve(length);
    
    std::uniform_int_distribution<> dist(0, 3);
    for (int i = 0; i < length; ++i) {
        kmer += nucleotides[dist(m_rng)];
    }
    
    return kmer;
} 