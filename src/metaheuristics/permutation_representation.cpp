#include "../../include/metaheuristics/representation.h"
#include "../../include/utils/random.h"
#include <algorithm>
#include <numeric>
#include <unordered_set>

std::vector<std::shared_ptr<Individual>> PermutationRepresentation::initializePopulation(
    size_t populationSize, const DNAInstance& instance) {
    if (populationSize == 0) {
        LOG_ERROR("Invalid population size: population size must be positive");
        throw std::invalid_argument("Population size must be positive");
    }
    
    std::vector<std::shared_ptr<Individual>> population;
    population.reserve(populationSize);
    
    // Simple initialization loop
    while (population.size() < populationSize) {
        auto individual = std::make_shared<Individual>();
        if (initializeIndividual(*individual, instance)) {
            population.push_back(individual);
        }
    }
    
    LOG_INFO("Generated initial population with {} individuals", population.size());
    return population;
}

bool PermutationRepresentation::initializeIndividual(Individual& individual, const DNAInstance& instance) {
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("Cannot initialize individual: spectrum is empty");
        return false;
    }

    const size_t spectrumSize = spectrum.size();
    
    // Create and shuffle indices
    std::vector<int> indices(spectrumSize);
    std::iota(indices.begin(), indices.end(), 0);
    
    auto& rng = Random::instance();
    for (size_t i = spectrumSize; i > 1; --i) {
        size_t j = rng.getRandomSizeT(0, i - 1);
        std::swap(indices[i - 1], indices[j]);
    }

    individual.setGenes(indices);
    return true;  // Accept any valid permutation
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    
    if (!individual) {
        return false;
    }
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) {
        return false;
    }
    
    // Simple validation - just check indices and DNA generation
    return validateGenes(genes, instance);
}

bool PermutationRepresentation::validateGenes(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    
    if (genes.empty()) return false;
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) return false;
    
    // Only check that genes are valid indices
    for (const auto& gene : genes) {
        if (gene < 0 || static_cast<size_t>(gene) >= spectrum.size()) return false;
    }
    
    // Generate DNA sequence
    std::string dna = toDNA(std::make_shared<Individual>(genes), instance);
    if (dna.empty()) return false;
    
    return true;
}

std::string PermutationRepresentation::toString(
    const std::shared_ptr<Individual>& individual,
    [[maybe_unused]] const DNAInstance& instance) const {
    if (!individual) return "";
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return "";

    std::string result;
    result.reserve(genes.size() * 3);  // Estimate 3 chars per number
    
    for (size_t i = 0; i < genes.size(); ++i) {
        if (i > 0) result += " ";
        result += std::to_string(genes[i]);
    }
    
    return result;
}

// Helper function to count mismatches between two strings
size_t PermutationRepresentation::countMismatches(const std::string& str1, const std::string& str2, size_t overlapLen) const {
    if (str1.empty() || str2.empty() || overlapLen == 0 || 
        overlapLen > str1.length() || overlapLen > str2.length()) {
        return SIZE_MAX;  // Invalid overlap
    }
    
    std::string suffix = str1.substr(str1.length() - overlapLen);
    std::string prefix = str2.substr(0, overlapLen);
    
    size_t mismatches = 0;
    for (size_t i = 0; i < overlapLen; ++i) {
        if (suffix[i] != prefix[i]) {
            mismatches++;
        }
    }
    
    return mismatches;
}

std::string PermutationRepresentation::toDNA(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const {
    if (!individual) return "";
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return "";
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) return "";
    
    size_t k = instance.getK();
    size_t deltaK = instance.getDeltaK();
    
    // No hard-coded length restrictions
    size_t targetLength = spectrum.size() + k;
    size_t maxLength = targetLength * 2;  // Just a reasonable upper bound
    
    std::string assembled;
    assembled.reserve(maxLength);
    
    // Start with first k-mer
    assembled = spectrum[genes[0]];
    
    for (size_t i = 1; i < genes.size(); ++i) {
        if (genes[i] < 0 || static_cast<size_t>(genes[i]) >= spectrum.size()) {
            continue;
        }
        
        const std::string& current = spectrum[genes[i]];
        bool foundOverlap = false;
        
        // Try overlaps from k-1 down to 1
        for (size_t overlapLen = k-1; overlapLen >= 1; --overlapLen) {
            size_t mismatches = countMismatches(assembled, current, overlapLen);
            if (mismatches <= deltaK) {
                assembled += current.substr(overlapLen);
                foundOverlap = true;
                break;
            }
        }
        
        // If no overlap found, just append with minimal connection
        if (!foundOverlap) {
            assembled += current.substr(1);
        }
    }
    
    return assembled;
} 