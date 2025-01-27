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
    
    // Keep track of unique solutions using their gene sequences
    std::vector<std::vector<int>> existingGenes;
    existingGenes.reserve(populationSize);
    
    int maxAttempts = populationSize * 5;  // Increased attempts to ensure success
    int attempts = 0;
    int consecutiveFailures = 0;
    const int MAX_CONSECUTIVE_FAILURES = 10;
    
    while (population.size() < static_cast<size_t>(populationSize) && attempts < maxAttempts) {
        auto individual = std::make_shared<Individual>();
        if (initializeIndividual(*individual, instance)) {
            const auto& genes = individual->getGenes();
            
            // Check if this is a duplicate solution
            bool isDuplicate = false;
            for (const auto& existing : existingGenes) {
                if (genes.size() == existing.size()) {
                    // Calculate similarity (Hamming distance)
                    size_t differences = 0;
                    for (size_t i = 0; i < genes.size(); ++i) {
                        if (genes[i] != existing[i]) {
                            differences++;
                        }
                    }
                    // If solutions are too similar (less than 20% different), consider it a duplicate
                    if (differences < genes.size() * 0.2) {
                        isDuplicate = true;
                        break;
                    }
                }
            }
            
            if (!isDuplicate) {
                population.push_back(individual);
                existingGenes.push_back(genes);
                consecutiveFailures = 0;  // Reset failure counter on success
            } else {
                consecutiveFailures++;
            }
        } else {
            consecutiveFailures++;
        }
        
        attempts++;
        
        // If too many consecutive failures, try to copy and modify existing solutions
        if (consecutiveFailures >= MAX_CONSECUTIVE_FAILURES && !population.empty()) {
            LOG_WARNING("Too many consecutive failures, using existing solutions as templates");
            
            // Take a random existing solution and modify it
            auto& rng = Random::instance();
            size_t idx = rng.getRandomInt(0, static_cast<int>(population.size() - 1));
            auto baseIndividual = population[idx];
            auto newIndividual = std::make_shared<Individual>(*baseIndividual);
            auto& genes = newIndividual->getGenes();
            
            // Perform more aggressive modifications
            size_t numSwaps = std::max(3UL, genes.size() / 3);  // Swap at least 33% of genes
            for (size_t i = 0; i < numSwaps; ++i) {
                size_t pos1 = rng.getRandomInt(0, static_cast<int>(genes.size() - 1));
                size_t pos2 = rng.getRandomInt(0, static_cast<int>(genes.size() - 1));
                std::swap(genes[pos1], genes[pos2]);
            }
            
            // Validate the modified individual
            if (validateGenes(genes, instance)) {
                population.push_back(newIndividual);
                existingGenes.push_back(genes);
                consecutiveFailures = 0;
            }
        }
    }
    
    // If we still don't have enough individuals, fill remaining slots with modified copies
    if (population.empty()) {
        LOG_ERROR("Failed to generate any valid individuals after " + std::to_string(attempts) + " attempts");
        throw std::runtime_error("Could not initialize population");
    }
    
    while (population.size() < static_cast<size_t>(populationSize)) {
        // Take a random existing solution and modify it more aggressively
        auto& rng = Random::instance();
        size_t idx = rng.getRandomInt(0, static_cast<int>(population.size() - 1));
        auto newIndividual = std::make_shared<Individual>(*population[idx]);
        auto& genes = newIndividual->getGenes();
        
        // Perform very aggressive modifications
        size_t numSwaps = std::max(4UL, genes.size() / 2);  // Swap at least 50% of genes
        for (size_t i = 0; i < numSwaps; ++i) {
            size_t pos1 = rng.getRandomInt(0, static_cast<int>(genes.size() - 1));
            size_t pos2 = rng.getRandomInt(0, static_cast<int>(genes.size() - 1));
            std::swap(genes[pos1], genes[pos2]);
        }
        
        population.push_back(newIndividual);
    }
    
    LOG_INFO("Generated initial population with " + std::to_string(population.size()) + 
             " individuals after " + std::to_string(attempts) + " attempts");
    return population;
}

bool PermutationRepresentation::initializeIndividual(Individual& individual, const DNAInstance& instance) {
    auto& rng = Random::instance();
    const auto& spectrum = instance.getSpectrum();
    
    if (spectrum.empty()) {
        LOG_ERROR("Empty spectrum");
        return false;
    }
    
    // Use spectrum size as the fixed length for all individuals
    size_t size = spectrum.size();
    
    // Create a permutation of indices
    std::vector<int> indices(size);
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ...
    
    // Shuffle the indices
    for (size_t i = indices.size() - 1; i > 0; --i) {
        size_t j = rng.getRandomInt(0, static_cast<int>(i));
        std::swap(indices[i], indices[j]);
    }
    
    // Keep track of used indices to avoid duplicates
    std::vector<bool> usedIndices(size, false);
    std::vector<int> finalIndices;
    finalIndices.reserve(size);  // Pre-allocate to avoid reallocation
    
    // Add all shuffled indices
    for (size_t i = 0; i < size; i++) {
        if (indices[i] >= 0 && static_cast<size_t>(indices[i]) < spectrum.size()) {
            finalIndices.push_back(indices[i]);
            usedIndices[indices[i]] = true;
        }
    }
    
    // Final validation
    if (finalIndices.empty() || finalIndices.size() != size) {
        LOG_ERROR("Failed to generate valid individual: size=" + std::to_string(finalIndices.size()) +
                  ", required size=" + std::to_string(size));
        return false;
    }
    
    // Validate all indices are within bounds
    for (int idx : finalIndices) {
        if (idx < 0 || static_cast<size_t>(idx) >= spectrum.size()) {
            LOG_ERROR("Invalid index generated: " + std::to_string(idx) + 
                     " (spectrum size: " + std::to_string(spectrum.size()) + ")");
            return false;
        }
    }
    
    individual.setGenes(finalIndices);
    individual.validate();  // Explicitly validate the individual
    return true;
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    
    if (!individual) return false;
    const auto& genes = individual->getGenes();
    if (genes.empty()) return false;
    
    // Use the more flexible validation
    return validateGenes(genes, instance);
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
    if (spectrum.empty()) return "";
    
    // Check if first gene is valid
    if (genes[0] < 0 || static_cast<size_t>(genes[0]) >= spectrum.size()) {
        LOG_ERROR("Invalid first gene index: " + std::to_string(genes[0]));
        return "";
    }
    
    const int k = instance.getK();
    std::string dna = spectrum[genes[0]];  // Start with first k-mer
    
    // Track used k-mers to ensure we're using the spectrum correctly
    std::vector<bool> used(spectrum.size(), false);
    used[genes[0]] = true;
    
    // For each remaining gene
    for (size_t i = 1; i < genes.size(); i++) {
        // Validate gene index
        if (genes[i] < 0 || static_cast<size_t>(genes[i]) >= spectrum.size()) {
            LOG_ERROR("Invalid gene index at position " + std::to_string(i) + ": " + std::to_string(genes[i]));
            return "";  // Invalid sequence
        }
        
        // Skip if already used and repetitions not allowed
        if (!instance.isRepAllowed() && used[genes[i]]) {
            LOG_DEBUG("Skipping repeated k-mer at position " + std::to_string(i));
            continue;
        }
        
        const std::string& current = spectrum[genes[i]];
        
        // Try different overlap lengths
        int bestOverlap = -1;
        int minMismatches = k + 1;  // More than maximum possible
        
        // Try overlaps from k-1 down to k/2
        for (int overlap = k - 1; overlap >= k/2; overlap--) {
            if (static_cast<int>(dna.length()) < overlap) continue;
            
            // Count mismatches in overlap region
            int mismatches = 0;
            for (int j = 0; j < overlap; j++) {
                if (dna[dna.length() - overlap + j] != current[j]) {
                    mismatches++;
                }
            }
            
            // Only accept overlaps with at most deltaK mismatches
            if (mismatches <= instance.getDeltaK() && mismatches < minMismatches) {
                minMismatches = mismatches;
                bestOverlap = overlap;
                
                // If we found a perfect match, no need to check smaller overlaps
                if (mismatches == 0) break;
            }
        }
        
        // Only append if we found a valid overlap
        if (bestOverlap > 0 && minMismatches <= instance.getDeltaK()) {
            dna += current.substr(bestOverlap);
            used[genes[i]] = true;
        } else {
            LOG_DEBUG("No valid overlap found for k-mer at position " + std::to_string(i));
            return "";  // Return empty string if no valid overlap found
        }
    }
    
    // Verify final sequence length
    if (dna.length() < instance.getK()) {
        LOG_DEBUG("Final sequence too short: " + std::to_string(dna.length()));
        return "";
    }
    
    // Verify that the sequence contains enough k-mers from the spectrum
    size_t usedCount = std::count(used.begin(), used.end(), true);
    if (usedCount < spectrum.size() * 0.7) {  // Require at least 70% coverage
        LOG_DEBUG("Insufficient spectrum coverage: " + std::to_string(usedCount) + "/" + std::to_string(spectrum.size()));
        return "";
    }
    
    return dna;
}

bool PermutationRepresentation::validateGenes(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    if (genes.empty()) return false;
    
    const auto& spectrum = instance.getSpectrum();
    
    // First check: all genes must be within valid bounds
    for (int gene : genes) {
        if (gene < 0 || static_cast<size_t>(gene) >= spectrum.size()) {
            return false;
        }
    }
    
    // Second check: genes must be the same size as spectrum
    // since we initialize with full spectrum size
    if (genes.size() != spectrum.size()) {
        return false;
    }
    
    return true;
} 