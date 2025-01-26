#include "../../include/dna/spectrum_generator.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

SpectrumGenerator::SpectrumGenerator(std::shared_ptr<std::mt19937> random)
    : m_random(random) {
    if (!random) {
        throw std::invalid_argument("Random generator cannot be null");
    }
}

std::vector<std::string> SpectrumGenerator::generateSpectrum(const std::string& dna, int k, int deltaK) {
    if (dna.empty()) {
        throw std::invalid_argument("DNA sequence cannot be empty");
    }
    if (k <= 0 || k > static_cast<int>(dna.length())) {
        throw std::invalid_argument("k-mer length must be positive and not greater than DNA length");
    }
    if (deltaK < 0) {
        throw std::invalid_argument("Delta K must be non-negative");
    }

    std::unordered_set<std::string> uniqueKmers;
    std::vector<std::string> kmers;

    // First k-mer always has length k
    uniqueKmers.insert(dna.substr(0, k));

    // Generate k-mers with variable length
    for (size_t i = 1; i <= dna.length() - k; ++i) {
        int currentK = k;
        
        // For all positions except the last k+2, allow length variation
        if (i < dna.length() - k - 2) {
            if (deltaK > 0) {
                int variation = Random::instance().getRandomInt(0, deltaK);
                if (variation > 0) {
                    // 50% chance for longer, 50% for shorter
                    if (Random::instance().generateProbability() < 0.5) {
                        currentK += variation;
                    } else {
                        currentK -= variation;
                    }
                }
            }
        }

        // Ensure we don't exceed DNA length
        if (i + currentK <= dna.length()) {
            std::string kmer = dna.substr(i, currentK);
            if (uniqueKmers.insert(kmer).second) {
                kmers.push_back(kmer);
            }
        }
    }

    // Sort k-mers lexicographically
    std::sort(kmers.begin(), kmers.end());
    return kmers;
}

void SpectrumGenerator::generateErrorKmers(const std::string& kmer, int errorsLeft, std::vector<std::string>& result) {
    if (errorsLeft <= 0) {
        return;
    }

    const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::uniform_int_distribution<int> posDist(0, kmer.length() - 1);
    std::uniform_int_distribution<int> nucDist(0, 3);

    // Generate a k-mer with one substitution error
    std::string errorKmer = kmer;
    int pos = posDist(*m_random);
    char original = errorKmer[pos];
    
    do {
        errorKmer[pos] = nucleotides[nucDist(*m_random)];
    } while (errorKmer[pos] == original);

    result.push_back(errorKmer);

    // Recursively generate k-mers with remaining errors
    if (errorsLeft > 1) {
        generateErrorKmers(errorKmer, errorsLeft - 1, result);
    }
} 