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
    : m_lNeg(lNeg)
    , m_random(std::random_device{}())
{
}

void NegativeErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (!validateInstance(instance)) {
        throw std::runtime_error("Cannot introduce negative errors: invalid instance");
    }

    if (instance.getLNeg() <= 0) {
        LOG_WARNING("No negative errors to introduce (lNeg <= 0)");
        return;
    }

    // Get a copy of the current spectrum
    auto spectrum = instance.getSpectrum();
    const size_t originalSize = spectrum.size();
    
    if (spectrum.empty()) {
        LOG_ERROR("Cannot introduce negative errors: spectrum is empty");
        throw std::runtime_error("Cannot introduce negative errors: spectrum is empty");
    }
    
    // Check if removing lNeg elements would make spectrum empty
    if (instance.getLNeg() >= static_cast<int>(spectrum.size())) {
        LOG_ERROR("Cannot remove more elements than spectrum size");
        throw std::runtime_error("Cannot remove more elements than spectrum size");
    }
    
    LOG_INFO("Introducing " + std::to_string(instance.getLNeg()) + " negative errors");
    LOG_INFO("Original spectrum size: " + std::to_string(originalSize));
    
    // Remove random k-mers from the spectrum
    std::uniform_int_distribution<> dist(0, spectrum.size() - 1);
    
    for (int i = 0; i < instance.getLNeg(); ++i) {
        int index = dist(m_random);
        spectrum.erase(spectrum.begin() + index);
        // Update distribution range after removing an element
        dist = std::uniform_int_distribution<>(0, spectrum.size() - 1);
    }
    
    // Update the instance's spectrum in a thread-safe way
    instance.setSpectrum(spectrum);
    LOG_INFO("Final spectrum size after negative errors: " + std::to_string(spectrum.size()));
}

PositiveErrorIntroducer::PositiveErrorIntroducer(int lPoz, int k)
    : m_lPoz(lPoz)
    , m_k(k)
    , m_random(std::random_device{}())
{
}

std::string PositiveErrorIntroducer::generateRandomKmer(int length) {
    const std::string nucleotides = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    std::string kmer;
    kmer.reserve(length);
    
    for (int i = 0; i < length; ++i) {
        kmer += nucleotides[dist(m_random)];
    }
    
    return kmer;
}

void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (!validateInstance(instance)) {
        throw std::runtime_error("Cannot introduce positive errors: invalid instance");
    }

    if (instance.getLPoz() <= 0) {
        LOG_WARNING("No positive errors to introduce (lPoz <= 0)");
        return;
    }

    // Get a copy of the current spectrum
    auto spectrum = instance.getSpectrum();
    const size_t originalSize = spectrum.size();
    
    if (spectrum.empty()) {
        LOG_ERROR("Cannot introduce positive errors: spectrum is empty");
        throw std::runtime_error("Cannot introduce positive errors: spectrum is empty");
    }
    
    LOG_INFO("Introducing " + std::to_string(instance.getLPoz()) + " positive errors");
    LOG_INFO("Original spectrum size: " + std::to_string(originalSize));
    
    // Generate random k-mers that are not in the spectrum
    const std::string nucleotides = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    
    for (int i = 0; i < instance.getLPoz(); ++i) {
        std::string kmer;
        bool unique;
        do {
            kmer.clear();
            for (int j = 0; j < instance.getK(); ++j) {
                kmer += nucleotides[dist(m_random)];
            }
            unique = std::find(spectrum.begin(), spectrum.end(), kmer) == spectrum.end();
        } while (!unique);
        
        spectrum.push_back(kmer);
    }
    
    // Update the instance's spectrum in a thread-safe way
    if (!spectrum.empty()) {
        instance.setSpectrum(spectrum);
        LOG_INFO("Final spectrum size after positive errors: " + std::to_string(spectrum.size()));
    } else {
        LOG_ERROR("Cannot set empty spectrum after positive errors");
        throw std::runtime_error("Cannot set empty spectrum after positive errors");
    }
}

std::unique_ptr<IErrorIntroductionStrategy> ErrorIntroducerFactory::createNegativeErrorIntroducer(int lNeg) {
    return std::make_unique<NegativeErrorIntroducer>(lNeg);
}

std::unique_ptr<IErrorIntroductionStrategy> ErrorIntroducerFactory::createPositiveErrorIntroducer(int lPoz, int k) {
    return std::make_unique<PositiveErrorIntroducer>(lPoz, k);
} 