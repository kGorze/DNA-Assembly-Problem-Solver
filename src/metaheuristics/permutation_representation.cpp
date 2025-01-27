#include "../../include/metaheuristics/representation.h"
#include "../../include/utils/random.h"
#include <algorithm>
#include <numeric>

std::vector<std::shared_ptr<Individual>> PermutationRepresentation::initializePopulation(
    int populationSize, const DNAInstance& instance) {
    if (populationSize <= 0) {
        throw std::invalid_argument("Population size must be positive");
    }
    
    std::vector<std::shared_ptr<Individual>> population;
    population.reserve(populationSize);
    
    for (int i = 0; i < populationSize; ++i) {
        auto individual = std::make_shared<Individual>();
        if (initializeIndividual(*individual, instance)) {
            population.push_back(individual);
        }
    }
    
    return population;
}

bool PermutationRepresentation::initializeIndividual(Individual& individual, const DNAInstance& instance) {
    auto& rng = Random::instance();
    const auto& spectrum = instance.getSpectrum();
    
    if (spectrum.empty()) {
        LOG_ERROR("Empty spectrum");
        return false;
    }
    
    // Create a simple permutation of indices
    std::vector<int> indices(spectrum.size());
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ...
    
    // Shuffle the indices
    for (size_t i = indices.size() - 1; i > 0; --i) {
        size_t j = rng.getRandomInt(0, static_cast<int>(i));
        std::swap(indices[i], indices[j]);
    }
    
    // Take a random subset of the indices (at least 2 elements)
    size_t size = rng.getRandomInt(2, static_cast<int>(indices.size()));
    indices.resize(size);
    
    individual.setGenes(indices);
    return true;  // Accept all valid permutations
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    
    if (!individual) return false;
    const auto& genes = individual->getGenes();
    if (genes.empty()) return false;
    
    // Only check if indices are within bounds
    for (int gene : genes) {
        if (gene < 0 || gene >= static_cast<int>(instance.getSpectrum().size())) {
            return false;
        }
    }
    
    return true;  // Accept any solution with valid indices
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

std::string PermutationRepresentation::toDNA(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    if (!individual) return "";
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return "";
    
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    std::string dna = spectrum[genes[0]];  // Start with first k-mer
    
    // Track used k-mers to avoid duplicates
    std::vector<bool> used(spectrum.size(), false);
    used[genes[0]] = true;
    
    // For each remaining gene
    for (size_t i = 1; i < genes.size(); i++) {
        if (static_cast<size_t>(genes[i]) >= spectrum.size()) continue;  // Skip invalid indices
        
        // Skip if already used and repetitions not allowed
        if (!instance.isRepAllowed() && used[genes[i]]) continue;
        
        const std::string& current = spectrum[genes[i]];
        
        // Try different overlap lengths
        int bestOverlap = -1;
        int minMismatches = k + 1;  // More than maximum possible
        
        // Try overlaps from k-1 down to 1
        for (int overlap = k - 1; overlap > 0; overlap--) {
            if (static_cast<int>(dna.length()) < overlap) continue;
            
            // Count mismatches in overlap region
            int mismatches = 0;
            for (int j = 0; j < overlap; j++) {
                if (dna[dna.length() - overlap + j] != current[j]) {
                    mismatches++;
                }
            }
            
            // If this is the best overlap so far, store it
            if (mismatches < minMismatches) {
                minMismatches = mismatches;
                bestOverlap = overlap;
            }
            
            // If we found a perfect match, no need to check smaller overlaps
            if (mismatches == 0) break;
        }
        
        // If we found any valid overlap, append with that overlap
        if (bestOverlap > 0) {
            // Even if there are mismatches, we'll still append but track the quality
            dna += current.substr(bestOverlap);
            used[genes[i]] = true;
        } else {
            // No overlap found - append with a single base overlap
            // This ensures we keep building the solution even with poor overlaps
            dna += current.substr(1);
            used[genes[i]] = true;
        }
    }
    
    return dna;  // Return the DNA sequence as is, without length adjustments
}

bool PermutationRepresentation::validateGenes(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    if (genes.empty()) return false;
    
    const auto& spectrum = instance.getSpectrum();
    
    // Only check if indices are within bounds
    for (int gene : genes) {
        if (gene < 0 || gene >= static_cast<int>(spectrum.size())) {
            return false;
        }
    }
    
    return true;  // Accept any solution with valid indices
} 