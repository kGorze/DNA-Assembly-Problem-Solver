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
    
    // Calculate minimum and maximum lengths based on instance parameters
    size_t minLength = std::max(
        static_cast<size_t>(spectrum.size() * 0.7),  // At least 70% of spectrum size
        static_cast<size_t>(instance.getK() + instance.getLPoz())  // Or k + lPoz
    );
    size_t maxLength = std::min(
        spectrum.size(),  // Full spectrum size
        minLength + static_cast<size_t>(instance.getLNeg())  // Or min length + lNeg for noise tolerance
    );
    
    // Check if valid length range exists
    if (minLength > maxLength || minLength > spectrum.size()) {
        LOG_ERROR("Invalid length range: min=" + std::to_string(minLength) + 
                  ", max=" + std::to_string(maxLength) + 
                  ", spectrum size=" + std::to_string(spectrum.size()));
        return false;
    }
    
    // Create a permutation of indices
    std::vector<int> indices(spectrum.size());
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ...
    
    // Shuffle the indices
    for (size_t i = indices.size() - 1; i > 0; --i) {
        size_t j = rng.getRandomInt(0, static_cast<int>(i));
        std::swap(indices[i], indices[j]);
    }
    
    // Generate a length that ensures good coverage
    size_t size = rng.getRandomInt(
        static_cast<int>(minLength),
        static_cast<int>(maxLength)
    );
    
    // Keep track of used indices to avoid duplicates
    std::vector<bool> usedIndices(spectrum.size(), false);
    std::vector<int> finalIndices;
    finalIndices.reserve(maxLength);
    
    // First add the initial shuffled sequence up to size
    for (size_t i = 0; i < size && i < indices.size(); i++) {
        finalIndices.push_back(indices[i]);
        usedIndices[indices[i]] = true;
    }
    
    // If we have noise parameters, ensure variety in solution length
    if (instance.getLNeg() > 0 || instance.getLPoz() > 0) {
        // 20% chance to add extra elements for handling negative errors
        if (rng.generateProbability() < 0.2 && finalIndices.size() < spectrum.size()) {
            size_t extraElements = rng.getRandomInt(1, std::min(instance.getLNeg(), 
                static_cast<int>(spectrum.size() - finalIndices.size())));
            
            // Try to add unused indices first
            std::vector<int> unusedIndices;
            for (size_t i = 0; i < spectrum.size(); i++) {
                if (!usedIndices[i]) {
                    unusedIndices.push_back(static_cast<int>(i));
                }
            }
            
            // Shuffle unused indices
            for (size_t i = unusedIndices.size() - 1; i > 0; --i) {
                size_t j = rng.getRandomInt(0, static_cast<int>(i));
                std::swap(unusedIndices[i], unusedIndices[j]);
            }
            
            // Add extra elements from unused indices
            for (size_t i = 0; i < extraElements && i < unusedIndices.size(); i++) {
                finalIndices.push_back(unusedIndices[i]);
                usedIndices[unusedIndices[i]] = true;
            }
        }
        // 20% chance to remove some elements for handling positive errors
        else if (rng.generateProbability() < 0.2 && finalIndices.size() > minLength) {
            size_t elementsToRemove = rng.getRandomInt(1, std::min(instance.getLPoz(), 
                static_cast<int>(finalIndices.size() - minLength)));
            finalIndices.resize(finalIndices.size() - elementsToRemove);
        }
    }
    
    // Ensure we have at least minLength elements
    if (finalIndices.size() < minLength) {
        // Add unused indices until we reach minLength
        for (size_t i = 0; i < spectrum.size() && finalIndices.size() < minLength; i++) {
            if (!usedIndices[i]) {
                finalIndices.push_back(i);
                usedIndices[i] = true;
            }
        }
        
        // If still not enough, allow reuse of indices
        while (finalIndices.size() < minLength) {
            int randomIndex = rng.getRandomInt(0, static_cast<int>(spectrum.size() - 1));
            finalIndices.push_back(randomIndex);
        }
    }
    
    // Final validation
    if (finalIndices.empty() || finalIndices.size() < minLength || finalIndices.size() > maxLength) {
        LOG_ERROR("Failed to generate valid individual: size=" + std::to_string(finalIndices.size()) +
                  ", required min=" + std::to_string(minLength) +
                  ", max=" + std::to_string(maxLength));
        return false;
    }
    
    individual.setGenes(finalIndices);
    return validateGenes(finalIndices, instance);
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
    
    // Only check if indices are within bounds
    for (int gene : genes) {
        if (gene < 0 || gene >= static_cast<int>(spectrum.size())) {
            return false;
        }
    }
    
    return true;  // Accept any solution with valid indices
} 