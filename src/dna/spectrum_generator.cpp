#include "../../include/dna/spectrum_generator.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

SpectrumGenerator::SpectrumGenerator(std::shared_ptr<std::mt19937> random)
    : m_random(random ? random : std::make_shared<std::mt19937>(std::random_device{}())) {}

std::vector<std::string> SpectrumGenerator::generateSpectrum(const std::string& dna, int k, int errorCount) {
    if (dna.empty()) {
        throw std::invalid_argument("DNA sequence cannot be empty");
    }
    if (k <= 0 || k > static_cast<int>(dna.length())) {
        throw std::invalid_argument("k-mer length must be positive and not greater than DNA length");
    }
    if (errorCount < 0) {
        throw std::invalid_argument("Error count cannot be negative");
    }

    std::unordered_set<std::string> uniqueKmers;

    // Ograniczamy liczbę błędów do jednego poziomu, niezależnie od przekazanego errorCount
    int effectiveErrorCount = (errorCount > 0) ? 1 : 0;

    // Generate all k-mers from the DNA sequence
    for (size_t i = 0; i <= dna.length() - k; ++i) {
        std::string kmer = dna.substr(i, k);
        uniqueKmers.insert(kmer);

        if (effectiveErrorCount > 0) {
            std::vector<std::string> errorKmers;
            generateErrorKmers(kmer, effectiveErrorCount, errorKmers);
            uniqueKmers.insert(errorKmers.begin(), errorKmers.end());
        }
    }

    // Convert set to vector and sort
    std::vector<std::string> result(uniqueKmers.begin(), uniqueKmers.end());
    std::sort(result.begin(), result.end());
    
    LOG_INFO("Generated spectrum with " + std::to_string(result.size()) + " k-mers");
    if (!result.empty()) {
        LOG_INFO("First k-mer: " + result.front() + ", Last k-mer: " + result.back());
    }
    
    return result;
}

void SpectrumGenerator::generateErrorKmers(const std::string& kmer, int errorsLeft, std::vector<std::string>& result) {
    // Zamiast rekurencyjnego generowania, generujemy tylko k-mery z jedną zamianą
    if (errorsLeft <= 0) {
        return;
    }

    const char nucleotides[] = {'A', 'C', 'G', 'T'};
    for (size_t pos = 0; pos < kmer.length(); ++pos) {
        std::string errorKmer = kmer;
        char original = errorKmer[pos];
        for (char nucleotide : nucleotides) {
            if (nucleotide != original) {
                errorKmer[pos] = nucleotide;
                result.push_back(errorKmer);
            }
        }
    }
} 