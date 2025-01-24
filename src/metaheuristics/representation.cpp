//
// Created by konrad_guest on 07/01/2025.
// SMART

#include "metaheuristics/representation.h"
#include "utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>

/************************************************
 *  DirectDNARepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
DirectDNARepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    int n = instance.getN();
    if (n <= 0) n = fallbackN; 

    std::mt19937 rng(std::random_device{}());
    static const char letters[] = "ACGT";

    for (int i = 0; i < popSize; i++) {
        auto individual = std::make_shared<std::vector<int>>();
        individual->reserve(n);

        for (int j = 0; j < n; j++) {
            int base = rng() % 4; // 0..3
            individual->push_back(base);
        }
        population.push_back(individual);
    }

    return population;
}

std::string 
DirectDNARepresentation::decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                                     const DNAInstance &instance)
{
    static const char letters[] = {'A','C','G','T'};
    if (!individual) return "";

    std::string dna;
    dna.reserve(individual->size());
    for (auto val : *individual) {
        int idx = std::max(0, std::min(val, 3));
        dna.push_back(letters[idx]);
    }

    return dna;
}

/************************************************
 *  PermutationRepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
PermutationRepresentation::initializePopulation(
    int populationSize,
    const DNAInstance& instance) 
{
    const size_t spectrumSize = instance.getSpectrum().size();
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        auto individual = std::make_shared<std::vector<int>>(spectrumSize);
        // Initialize with sequential numbers
        std::iota(individual->begin(), individual->end(), 0);
        // Shuffle the sequence
        std::random_shuffle(individual->begin(), individual->end());
        
        // Validate the individual
        bool valid = true;
        for (size_t j = 0; j < individual->size(); ++j) {
            if ((*individual)[j] < 0 || (*individual)[j] >= static_cast<int>(spectrumSize)) {
                valid = false;
                break;
            }
        }
        
        if (valid) {
            population.push_back(individual);
        } else {
            --i; // Try again for this position
        }
    }
    return population;
}

std::string PermutationRepresentation::decodeToDNA(
    std::shared_ptr<std::vector<int>> individual,
    const DNAInstance &instance)
{
    if (!individual || individual->empty()) {
        std::cerr << "[ERROR decodeToDNA] Individual is nullptr or empty\n";
        return "";
    }
    const auto &spec = instance.getSpectrum();
    int k = instance.getK();
    int n = instance.getN();
    if (k <= 0 || n <= 0) {
        std::cerr << "[ERROR decodeToDNA] Invalid k or n (k=" << k << ", n=" << n << ")\n";
        return "";
    }

    int firstIndex = (*individual)[0];
    if (firstIndex < 0 || firstIndex >= (int)spec.size()) {
        std::cerr << "[ERROR decodeToDNA] firstIndex out-of-range: " << firstIndex << "\n";
        return "";
    }

    std::string result = spec[firstIndex]; // ma długość co najmniej k (z założenia)

    for (size_t i = 1; i < individual->size(); i++) {
        int fragIdx = (*individual)[i];
        if (fragIdx < 0 || fragIdx >= (int)spec.size()) {
            std::cerr << "[ERROR decodeToDNA] fragIdx out-of-range: " << fragIdx << " at position " << i << "\n";
            result.push_back('N'); 
            continue;
        }

        const std::string &frag = spec[fragIdx];
        if ((int)frag.size() < k - 2) { // Allow for deltaK=2
            std::cerr << "[WARNING decodeToDNA] Fragment size " << frag.size() 
                     << " smaller than k-deltaK (" << k-2 << ")\n";
            continue;
        }

        // For variable length fragments, always take the last character
        result.push_back(frag.back());
    }

    // Ensure result length matches target length
    if ((int)result.size() > n) {
        std::cerr << "[WARNING decodeToDNA] Reconstructed DNA longer than n, resizing from " << result.size() << " to " << n << "\n";
        result.resize(n);
    } else while ((int)result.size() < n) {
        result.push_back('N');
    }

    // **Finalny Debug Log**: Wypisz pierwszy i ostatni kilka znaków DNA
    std::cout << "[DEBUG decodeToDNA] Reconstructed DNA (first 10 chars): " 
              << result.substr(0, 10) << "...\n";
    std::cout << "[DEBUG decodeToDNA] Reconstructed DNA (last 10 chars): " 
              << result.substr(result.size() - 10, 10) << "\n";

    return result;
}

/************************************************
 *  GraphPathRepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
GraphPathRepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    // Calculate required length based on n and k
    int k = instance.getK();
    int n = instance.getN();
    int requiredLength = n - k + 1;  // This ensures we get exactly n characters after decoding
    
    for(int i = 0; i < popSize; i++) {
        auto path = std::make_shared<std::vector<int>>();
        path->reserve(requiredLength);
        
        // Fill with indices from 0 to spectrum.size()-1
        for(int idx = 0; idx < (int)instance.getSpectrum().size(); idx++) {
            path->push_back(idx);
        }
        
        // Shuffle and resize to required length
        std::shuffle(path->begin(), path->end(), std::mt19937(std::random_device{}()));
        path->resize(requiredLength);
        
        population.push_back(path);
    }

    return population;
}

std::string 
GraphPathRepresentation::decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                                     const DNAInstance &instance)
{
    if (!individual || individual->empty()) {
        return "";
    }
    const auto &spec = instance.getSpectrum();
    if (spec.empty()) {
        return "";
    }

    int k = instance.getK();
    int n = instance.getN();
    
    // Get first fragment
    int firstIndex = (*individual)[0];
    if (firstIndex < 0 || firstIndex >= (int)spec.size()) {
        return "";
    }

    std::string result = spec[firstIndex];

    // Add remaining characters
    for (size_t i = 1; i < individual->size(); i++) {
        int fragIdx = (*individual)[i];
        if (fragIdx < 0 || fragIdx >= (int)spec.size()) {
            continue;
        }

        const std::string &frag = spec[fragIdx];
        if ((int)frag.size() < k - 2) { // Allow for deltaK=2
            std::cerr << "[WARNING decodeToDNA] Fragment size " << frag.size() 
                     << " smaller than k-deltaK (" << k-2 << ")\n";
            continue;
        }
        result.push_back(frag.back());
    }

    // Result should be exactly n characters
    if ((int)result.size() != n) {
        LOG_WARNING("Decoded DNA length (" + std::to_string(result.size()) + 
                   ") differs from required length (" + std::to_string(n) + ")");
        return "";  // Return empty string to indicate invalid solution
    }

    return result;
}
