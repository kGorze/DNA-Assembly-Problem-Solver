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
    
    // Always use the full spectrum size
    size_t size = spectrum.size();
    
    // Create a permutation of indices
    std::vector<int> indices(size);
    std::iota(indices.begin(), indices.end(), 0);  // Fill with 0, 1, 2, ...
    
    // Shuffle the indices
    for (size_t i = indices.size() - 1; i > 0; --i) {
        size_t j = rng.getRandomInt(0, static_cast<int>(i));
        std::swap(indices[i], indices[j]);
    }
    
    // Validate indices
    for (int index : indices) {
        if (index < 0 || static_cast<size_t>(index) >= spectrum.size()) {
            LOG_ERROR("Generated invalid index: " + std::to_string(index) + 
                      ", spectrum size: " + std::to_string(spectrum.size()));
            return false;
        }
    }
    
    // Set the genes
    individual.setGenes(indices);
    LOG_DEBUG("Initialized individual with " + std::to_string(indices.size()) + " genes");
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
            continue;  // Skip invalid genes instead of returning empty string
        }
        
        // Skip if already used and repetitions not allowed
        if (!instance.isRepAllowed() && used[genes[i]]) {
            LOG_DEBUG("Skipping repeated k-mer at position " + std::to_string(i));
            continue;
        }
        
        const std::string& current = spectrum[genes[i]];
        
        // Try different overlap lengths
        int bestOverlap = -1;
        int minMismatches = k;  // Allow up to k mismatches (more lenient)
        
        // Try overlaps from k-1 down to 1 (more lenient minimum overlap)
        for (int overlap = k - 1; overlap >= 1; overlap--) {
            if (static_cast<int>(dna.length()) < overlap) continue;
            
            // Count mismatches in overlap region
            int mismatches = 0;
            for (int j = 0; j < overlap; j++) {
                if (dna[dna.length() - overlap + j] != current[j]) {
                    mismatches++;
                }
            }
            
            // Accept overlaps with reasonable number of mismatches
            // More lenient: accept if mismatches <= min(deltaK * 2, overlap/2)
            int maxAllowedMismatches = std::min(instance.getDeltaK() * 2, overlap / 2);
            if (mismatches <= maxAllowedMismatches && mismatches < minMismatches) {
                minMismatches = mismatches;
                bestOverlap = overlap;
                
                // If we found a good enough match, no need to check smaller overlaps
                if (mismatches <= instance.getDeltaK()) break;
            }
        }
        
        // If no good overlap found, try minimal overlap
        if (bestOverlap <= 0) {
            // Use minimal overlap of 1 base
            bestOverlap = 1;
            dna += current.substr(bestOverlap);
        } else {
            dna += current.substr(bestOverlap);
        }
        used[genes[i]] = true;
    }
    
    // Verify final sequence length - more lenient minimum length
    if (dna.length() < instance.getK() / 2) {  // Allow sequences at least half the k-mer length
        LOG_DEBUG("Final sequence too short: " + std::to_string(dna.length()));
        return "";
    }
    
    // More lenient coverage requirement - accept if we used at least 30% of k-mers
    size_t usedCount = std::count(used.begin(), used.end(), true);
    if (usedCount < spectrum.size() * 0.3) {  // Reduced from 70% to 30%
        LOG_DEBUG("Low spectrum coverage: " + std::to_string(usedCount) + "/" + std::to_string(spectrum.size()));
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
            LOG_DEBUG("Invalid gene value: " + std::to_string(gene) + 
                     " (spectrum size: " + std::to_string(spectrum.size()) + ")");
            return false;
        }
    }
    
    // Second check: genes size should be reasonable
    // Allow sizes between 50% of spectrum size and full spectrum size
    size_t minSize = std::max(instance.getK(), static_cast<int>(spectrum.size() / 2));
    if (genes.size() < minSize || genes.size() > spectrum.size()) {
        LOG_DEBUG("Invalid genes size: " + std::to_string(genes.size()) + 
                 " (min: " + std::to_string(minSize) + 
                 ", max: " + std::to_string(spectrum.size()) + ")");
        return false;
    }
    
    // Third check: if repetitions are not allowed, check for duplicates
    if (!instance.isRepAllowed()) {
        std::vector<bool> used(spectrum.size(), false);
        for (int gene : genes) {
            if (used[gene]) {
                LOG_DEBUG("Duplicate gene found: " + std::to_string(gene));
                return false;
            }
            used[gene] = true;
        }
    }
    
    return true;
} 