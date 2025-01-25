#include "dna/error_introduction.h"
#include "utils/logging.h"
#include <algorithm>
#include <random>
#include <stdexcept>

bool BaseErrorIntroducer::validateInstance(const DNAInstance& instance) const {
    if (instance.getDNA().empty()) {
        LOG_ERROR("Cannot introduce errors: DNA sequence is empty");
        return false;
    }
    if (instance.getSpectrum().empty()) {
        LOG_ERROR("Cannot introduce errors: spectrum is empty");
        return false;
    }
    return true;
}

NegativeErrorIntroducer::NegativeErrorIntroducer(int lNeg)
    : m_lNeg(lNeg), m_random(std::random_device{}()) {}

void NegativeErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (instance.getLNeg() <= 0) {
        LOG_WARNING("No negative errors to introduce (lNeg <= 0)");
        return;
    }

    // Get current spectrum
    auto spectrum = instance.getSpectrum();
    
    // Generate random k-mers that are not in the spectrum
    const std::string nucleotides = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    
    for (int i = 0; i < instance.getLNeg(); ++i) {
        std::string kmer;
        bool unique;
        do {
            kmer.clear();
            for (int j = 0; j < instance.getK(); ++j) {
                kmer += nucleotides[dist(m_random->getGenerator())];
            }
            unique = std::find(spectrum.begin(), spectrum.end(), kmer) == spectrum.end();
        } while (!unique);
        
        spectrum.push_back(kmer);
    }
    
    // Update the instance's spectrum
    instance.setSpectrum(spectrum);
}

PositiveErrorIntroducer::PositiveErrorIntroducer(int lPoz, int k)
    : m_lPoz(lPoz), m_k(k), m_random(std::random_device{}()) {}

std::string PositiveErrorIntroducer::generateRandomKmer(int length) const {
    const std::string nucleotides = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    std::string kmer;
    kmer.reserve(length);
    
    for (int i = 0; i < length; ++i) {
        kmer += nucleotides[dist(m_random->getGenerator())];
    }
    
    return kmer;
}

void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (instance.getLPoz() <= 0) {
        LOG_WARNING("No positive errors to introduce (lPoz <= 0)");
        return;
    }

    // Get current spectrum
    auto spectrum = instance.getSpectrum();
    
    // Remove random k-mers from the spectrum
    std::uniform_int_distribution<> dist(0, spectrum.size() - 1);
    
    for (int i = 0; i < std::min(instance.getLPoz(), static_cast<int>(spectrum.size())); ++i) {
        int index = dist(m_random->getGenerator());
        spectrum.erase(spectrum.begin() + index);
    }
    
    // Update the instance's spectrum
    instance.setSpectrum(spectrum);
}

std::unique_ptr<IErrorIntroductionStrategy> ErrorIntroducerFactory::createNegativeErrorIntroducer(int lNeg) {
    return std::make_unique<NegativeErrorIntroducer>(lNeg);
}

std::unique_ptr<IErrorIntroductionStrategy> ErrorIntroducerFactory::createPositiveErrorIntroducer(int lPoz, int k) {
    return std::make_unique<PositiveErrorIntroducer>(lPoz, k);
} 