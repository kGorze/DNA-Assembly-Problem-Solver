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
        LOG_ERROR("Initialization failed: Empty spectrum");
        return false;
    }
    
    LOG_DEBUG("Initializing individual with spectrum size " + std::to_string(spectrum.size()));
    
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
    for (size_t i = 0; i < indices.size(); i++) {
        if (indices[i] < 0 || static_cast<size_t>(indices[i]) >= spectrum.size()) {
            LOG_ERROR("Generated invalid index at position " + std::to_string(i) + 
                     ": " + std::to_string(indices[i]) + 
                     " (spectrum size: " + std::to_string(spectrum.size()) + ")");
            return false;
        }
    }
    
    // Set the genes
    individual.setGenes(indices);
    LOG_DEBUG("Successfully initialized individual with " + std::to_string(indices.size()) + 
              " genes (spectrum size: " + std::to_string(spectrum.size()) + ")");
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
    // Basic validation remains the same
    if (instance.getK() <= 0) {
        LOG_ERROR("toDNA failed: Invalid k-mer length " + std::to_string(instance.getK()));
        return "";
    }
    if (instance.getDeltaK() < 0) {
        LOG_ERROR("toDNA failed: Invalid deltaK " + std::to_string(instance.getDeltaK()));
        return "";
    }
    if (!individual || individual->getGenes().empty()) {
        LOG_ERROR("toDNA failed: Invalid individual");
        return "";
    }

    const auto& genes = individual->getGenes();
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("toDNA failed: Empty spectrum");
        return "";
    }

    LOG_DEBUG("Starting DNA assembly with " + std::to_string(genes.size()) + " genes");
    
    // Initialize with first valid k-mer
    size_t startIdx = 0;
    while (startIdx < genes.size()) {
        if (genes[startIdx] >= 0 && static_cast<size_t>(genes[startIdx]) < spectrum.size()) {
            break;
        }
        startIdx++;
    }
    if (startIdx >= genes.size()) {
        LOG_ERROR("No valid starting k-mer found");
        return "";
    }

    std::string dna = spectrum[genes[startIdx]];
    const int k = instance.getK();
    std::vector<bool> used(spectrum.size(), false);
    used[genes[startIdx]] = true;

    // Dynamic overlap threshold based on k-mer length
    const int minOverlap = std::max(1, k / 4);  // More permissive minimum overlap
    const int maxMismatchesPerOverlap = std::max(instance.getDeltaK(), k / 3);  // More tolerant of mismatches

    // For each remaining gene
    for (size_t i = startIdx + 1; i < genes.size(); i++) {
        if (genes[i] < 0 || static_cast<size_t>(genes[i]) >= spectrum.size()) {
            LOG_DEBUG("Skipping invalid gene at position " + std::to_string(i));
            continue;
        }

        if (!instance.isRepAllowed() && used[genes[i]]) {
            LOG_DEBUG("Skipping repeated k-mer at position " + std::to_string(i));
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        int bestOverlap = -1;
        int bestMismatches = k;
        std::string bestExtension;

        // Try overlaps from largest to smallest
        for (int overlap = k - 1; overlap >= minOverlap; overlap--) {
            if (dna.length() < static_cast<size_t>(overlap)) continue;

            int mismatches = 0;
            for (int j = 0; j < overlap; j++) {
                if (dna[dna.length() - overlap + j] != current[j]) {
                    mismatches++;
                }
            }

            // Accept if mismatches are reasonable for this overlap length
            int maxAllowedMismatches = std::min(maxMismatchesPerOverlap, 
                                               static_cast<int>(overlap * 0.4));  // Allow up to 40% mismatches
            
            if (mismatches <= maxAllowedMismatches && 
                (bestOverlap == -1 || mismatches <= bestMismatches)) {
                bestOverlap = overlap;
                bestMismatches = mismatches;
                bestExtension = current.substr(overlap);
                
                // If we found a very good match, use it immediately
                if (mismatches <= instance.getDeltaK()) break;
            }
        }

        // If no good overlap found, try minimal overlap with relaxed criteria
        if (bestOverlap == -1) {
            LOG_DEBUG("Using minimal overlap for k-mer at position " + std::to_string(i));
            bestExtension = current.substr(minOverlap);
            bestOverlap = minOverlap;
        }

        dna += bestExtension;
        used[genes[i]] = true;
        
        LOG_DEBUG("Added k-mer with overlap " + std::to_string(bestOverlap) + 
                 ", mismatches " + std::to_string(bestMismatches) + 
                 ", new length " + std::to_string(dna.length()));
    }

    // More flexible minimum length requirement
    size_t minLength = std::max(static_cast<size_t>(k), spectrum.size() / 4);
    if (dna.length() < minLength) {
        LOG_DEBUG("Assembled sequence too short: " + std::to_string(dna.length()) + 
                 " < " + std::to_string(minLength));
        return "";
    }

    // Dynamic coverage requirement based on sequence length
    size_t usedCount = std::count(used.begin(), used.end(), true);
    double coveragePercent = (static_cast<double>(usedCount) / spectrum.size()) * 100;
    
    // Adjust minimum coverage based on sequence length
    double minCoverage = 20.0;  // Base minimum coverage
    if (dna.length() > spectrum.size() * k / 2) {
        minCoverage = 15.0;  // Lower requirement for longer sequences
    }
    
    LOG_DEBUG("Coverage: " + std::to_string(coveragePercent) + "% (minimum: " + 
              std::to_string(minCoverage) + "%)");
    
    if (coveragePercent < minCoverage) {
        LOG_DEBUG("Insufficient coverage");
        return "";
    }

    LOG_DEBUG("Successfully assembled DNA sequence of length " + std::to_string(dna.length()) + 
              " using " + std::to_string(usedCount) + " k-mers");
    return dna;
}

bool PermutationRepresentation::validateGenes(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    // Validate DNAInstance parameters first
    if (instance.getK() <= 0) {
        LOG_ERROR("Invalid k-mer length: " + std::to_string(instance.getK()) + " (must be positive)");
        return false;
    }
    if (instance.getDeltaK() < 0) {
        LOG_ERROR("Invalid deltaK: " + std::to_string(instance.getDeltaK()) + " (must be non-negative)");
        return false;
    }
    LOG_DEBUG("DNAInstance parameters: k=" + std::to_string(instance.getK()) + 
              ", deltaK=" + std::to_string(instance.getDeltaK()) +
              ", repetitions=" + std::string(instance.isRepAllowed() ? "allowed" : "not allowed"));

    if (genes.empty()) {
        LOG_DEBUG("Validation failed: Empty genes vector");
        return false;
    }
    
    const auto& spectrum = instance.getSpectrum();
    LOG_DEBUG("Validating genes vector of size " + std::to_string(genes.size()) + 
              " against spectrum size " + std::to_string(spectrum.size()));
    
    // First check: all genes must be within valid bounds
    for (size_t i = 0; i < genes.size(); i++) {
        if (genes[i] < 0 || static_cast<size_t>(genes[i]) >= spectrum.size()) {
            LOG_DEBUG("Validation failed: Invalid gene value at index " + std::to_string(i) + 
                     ": " + std::to_string(genes[i]) + 
                     " (valid range: [0, " + std::to_string(spectrum.size() - 1) + "])");
            return false;
        }
    }
    
    // Second check: genes size should be reasonable
    // Allow sizes between 50% of spectrum size and 150% of spectrum size
    size_t minSize = std::max(instance.getK(), static_cast<int>(spectrum.size() * 0.5));
    size_t maxSize = static_cast<size_t>(spectrum.size() * 1.5);
    LOG_DEBUG("Size validation: Current=" + std::to_string(genes.size()) + 
              ", Min=" + std::to_string(minSize) + 
              ", Max=" + std::to_string(maxSize));
    
    if (genes.size() < minSize || genes.size() > maxSize) {
        LOG_DEBUG("Validation failed: Invalid genes size " + std::to_string(genes.size()) + 
                 " (required range: [" + std::to_string(minSize) + ", " + 
                 std::to_string(maxSize) + "])");
        return false;
    }
    
    // Third check: if repetitions are not allowed, be more lenient with duplicates
    if (!instance.isRepAllowed()) {
        LOG_DEBUG("Checking for duplicates (repetitions not allowed)");
        // Count occurrences of each gene
        std::vector<int> geneCount(spectrum.size(), 0);
        for (size_t i = 0; i < genes.size(); i++) {
            geneCount[genes[i]]++;
            if (geneCount[genes[i]] > 2) {
                LOG_DEBUG("Validation failed: Gene " + std::to_string(genes[i]) + 
                         " appears more than twice (found at index " + std::to_string(i) + 
                         ", count=" + std::to_string(geneCount[genes[i]]) + ")");
                return false;
            }
        }
        
        // Log distribution of gene occurrences
        std::stringstream ss;
        ss << "Gene occurrence distribution: ";
        int singleCount = 0, doubleCount = 0;
        for (size_t i = 0; i < geneCount.size(); i++) {
            if (geneCount[i] == 1) singleCount++;
            else if (geneCount[i] == 2) doubleCount++;
        }
        ss << singleCount << " single, " << doubleCount << " double occurrences";
        LOG_DEBUG(ss.str());
    }
    
    LOG_DEBUG("Gene validation successful");
    return true;
} 