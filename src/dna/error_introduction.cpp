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
    return true;
}

void NegativeErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (!validateInstance(instance)) {
        throw std::runtime_error("Invalid instance for negative error introduction");
    }

    if (m_lNeg <= 0) {
        LOG_INFO("No negative errors to introduce");
        return;
    }

    auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("Cannot introduce negative errors: spectrum is empty");
        throw std::runtime_error("Empty spectrum");
    }

    std::uniform_int_distribution<size_t> dist(0, spectrum.size() - 1);
    std::vector<size_t> indicesToRemove;

    for (int i = 0; i < m_lNeg && !spectrum.empty(); ++i) {
        size_t index = dist(m_rng);
        indicesToRemove.push_back(index);
    }

    // Sort indices in descending order to remove from back to front
    std::sort(indicesToRemove.begin(), indicesToRemove.end(), std::greater<>());

    for (size_t index : indicesToRemove) {
        if (index < spectrum.size()) {
            spectrum.erase(spectrum.begin() + index);
        }
    }
}

std::string PositiveErrorIntroducer::generateRandomKmer(int length) const {
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<> dist(0, 3);
    std::string kmer;
    kmer.reserve(length);
    
    for (int i = 0; i < length; ++i) {
        kmer += nucleotides[dist(m_rng)];
    }
    
    return kmer;
}

void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (!validateInstance(instance)) {
        throw std::runtime_error("Invalid instance for positive error introduction");
    }

    if (m_lPoz <= 0) {
        LOG_INFO("No positive errors to introduce");
        return;
    }

    auto& spectrum = instance.getSpectrum();
    for (int i = 0; i < m_lPoz; ++i) {
        std::string randomKmer = generateRandomKmer(m_k);
        spectrum.push_back(randomKmer);
    }

    // Remove duplicates that might have been introduced
    std::sort(spectrum.begin(), spectrum.end());
    spectrum.erase(std::unique(spectrum.begin(), spectrum.end()), spectrum.end());
}

std::unique_ptr<IErrorIntroductionStrategy> ErrorIntroducerFactory::createNegativeErrorIntroducer(int lNeg) {
    return std::make_unique<NegativeErrorIntroducer>(lNeg);
}

std::unique_ptr<IErrorIntroductionStrategy> ErrorIntroducerFactory::createPositiveErrorIntroducer(int lPoz, int k) {
    return std::make_unique<PositiveErrorIntroducer>(lPoz, k);
} 