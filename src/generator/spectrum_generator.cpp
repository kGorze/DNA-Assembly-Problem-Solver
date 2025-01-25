#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <random>
#include <algorithm>

std::vector<std::string> SpectrumGenerator::generateSpectrum(const std::string& dna, int k, int deltaK) {
    if (dna.empty()) {
        LOG_ERROR("Cannot generate spectrum from empty DNA");
        throw std::invalid_argument("DNA sequence is empty");
    }
    if (k <= 0) {
        LOG_ERROR("Invalid k value: " + std::to_string(k));
        throw std::invalid_argument("k must be positive");
    }
    if (deltaK < 0) {
        LOG_ERROR("Invalid deltaK value: " + std::to_string(deltaK));
        throw std::invalid_argument("deltaK cannot be negative");
    }
    if (k > static_cast<int>(dna.length())) {
        LOG_ERROR("k value (" + std::to_string(k) + ") larger than DNA length (" + std::to_string(dna.length()) + ")");
        throw std::invalid_argument("k value larger than DNA length");
    }

    std::vector<std::string> spectrum;
    static thread_local std::random_device rd;
    static thread_local std::mt19937 gen(rd());
    std::uniform_int_distribution<> deltaKDist(0, deltaK);
    std::uniform_int_distribution<> signDist(0, 1);

    // First oligonucleotide always has length k
    spectrum.push_back(dna.substr(0, k));

    // Generate middle oligonucleotides with variable length
    size_t pos = 1;
    while (pos + k <= dna.length() - (k + 2)) {
        int delta = deltaKDist(gen);
        int length = k;
        if (delta > 0) {
            length += signDist(gen) ? delta : -delta;
        }
        
        // Ensure length is valid
        length = std::max(1, std::min(length, static_cast<int>(dna.length() - pos)));
        
        spectrum.push_back(dna.substr(pos, length));
        pos++;
    }

    // Last k+2 oligonucleotides have length k
    while (pos + k <= dna.length()) {
        spectrum.push_back(dna.substr(pos, k));
        pos++;
    }

    // Sort the spectrum
    std::sort(spectrum.begin(), spectrum.end());

    return spectrum;
} 