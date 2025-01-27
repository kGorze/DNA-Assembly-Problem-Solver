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
    
    std::vector<std::vector<size_t>> existingGenes;
    existingGenes.reserve(populationSize);
    
    // Calculate max attempts using size_t
    size_t maxAttempts = populationSize * 5;
    
    size_t attempts = 0;
    size_t consecutiveFailures = 0;
    const size_t MAX_CONSECUTIVE_FAILURES = 10;
    
    while (population.size() < populationSize && attempts < maxAttempts) {
        auto individual = std::make_shared<Individual>();
        if (initializeIndividual(*individual, instance)) {
            const auto& genes = individual->getGenes();
            
            // Compare differences using size_t
            bool isDuplicate = false;
            for (const auto& existing : existingGenes) {
                if (genes.size() == existing.size()) {
                    size_t differences = 0;
                    for (size_t i = 0; i < genes.size(); ++i) {
                        const size_t geneValue = static_cast<size_t>(genes[i]);
                        if (geneValue != existing[i]) {
                            differences++;
                        }
                    }
                    // Use size_t for threshold calculation
                    const size_t threshold = genes.size() / 5;  // 20% threshold
                    if (differences < threshold) {
                        isDuplicate = true;
                        break;
                    }
                }
            }
            
            if (!isDuplicate) {
                population.push_back(individual);
                // Convert genes to size_t before storing
                std::vector<size_t> sizeGenes;
                sizeGenes.reserve(genes.size());
                for (const int gene : genes) {
                    if (gene >= 0) {  // Only convert non-negative values
                        sizeGenes.push_back(static_cast<size_t>(gene));
                    }
                }
                existingGenes.push_back(sizeGenes);
                consecutiveFailures = 0;
            } else {
                consecutiveFailures++;
            }
        } else {
            consecutiveFailures++;
        }
        
        attempts++;
        
        if (consecutiveFailures >= MAX_CONSECUTIVE_FAILURES && !population.empty()) {
            // LOG_WARNING("Too many consecutive failures, using existing solutions as templates");  // Commented out - not essential
            
            auto& rng = Random::instance();
            if (population.empty()) {
                LOG_ERROR("Cannot generate random index: population is empty");
                continue;
            }
            
            const size_t maxIdx = population.size() - 1;
            const size_t idx = rng.getRandomSizeT(0, maxIdx);
            
            auto newIndividual = std::make_shared<Individual>(*population[idx]);
            auto& genes = newIndividual->getGenes();
            
            if (genes.empty()) {
                LOG_ERROR("Invalid individual: empty genes vector");
                continue;
            }
            
            const size_t maxGeneIdx = genes.size() - 1;
            const size_t numSwaps = std::max(size_t{3}, genes.size() / 3);
            
            for (size_t i = 0; i < numSwaps; ++i) {
                const size_t pos1 = rng.getRandomSizeT(0, maxGeneIdx);
                const size_t pos2 = rng.getRandomSizeT(0, maxGeneIdx);
                
                if (pos1 <= maxGeneIdx && pos2 <= maxGeneIdx) {
                    std::swap(genes[pos1], genes[pos2]);
                } else {
                    LOG_ERROR("Invalid swap positions generated: " + 
                             std::to_string(pos1) + ", " + std::to_string(pos2) + 
                             " (max: " + std::to_string(maxGeneIdx) + ")");
                    continue;
                }
            }
            
            if (validateGenes(genes, instance)) {
                population.push_back(newIndividual);
                // Convert genes to size_t before storing
                std::vector<size_t> sizeGenes;
                sizeGenes.reserve(genes.size());
                for (const int gene : genes) {
                    if (gene >= 0) {  // Only convert non-negative values
                        sizeGenes.push_back(static_cast<size_t>(gene));
                    }
                }
                existingGenes.push_back(sizeGenes);
                consecutiveFailures = 0;
            }
        }
    }
    
    // Fill remaining population slots
    while (population.size() < populationSize) {
        if (population.empty()) {
            LOG_ERROR("Cannot generate random index: population is empty");
            break;
        }
        
        auto& rng = Random::instance();
        const size_t maxIdx = population.size() - 1;
        const size_t idx = rng.getRandomSizeT(0, maxIdx);
        
        auto newIndividual = std::make_shared<Individual>(*population[idx]);
        auto& genes = newIndividual->getGenes();
        
        if (genes.empty()) {
            LOG_ERROR("Invalid individual: empty genes vector");
            continue;
        }
        
        const size_t maxGeneIdx = genes.size() - 1;
        const size_t numSwaps = std::max(size_t{4}, genes.size() / 2);
        
        for (size_t i = 0; i < numSwaps; ++i) {
            const size_t pos1 = rng.getRandomSizeT(0, maxGeneIdx);
            const size_t pos2 = rng.getRandomSizeT(0, maxGeneIdx);
            
            if (pos1 <= maxGeneIdx && pos2 <= maxGeneIdx) {
                std::swap(genes[pos1], genes[pos2]);
            } else {
                LOG_ERROR("Invalid swap positions generated: " + 
                         std::to_string(pos1) + ", " + std::to_string(pos2) + 
                         " (max: " + std::to_string(maxGeneIdx) + ")");
                continue;
            }
        }
        
        if (validateGenes(genes, instance)) {
            population.push_back(newIndividual);
        }
    }
    
    LOG_INFO("Generated initial population with " + std::to_string(population.size()) + 
             " individuals after " + std::to_string(attempts) + " attempts");
    return population;
}

bool PermutationRepresentation::initializeIndividual(Individual& individual, const DNAInstance& instance) {
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("Cannot initialize individual: spectrum is empty");
        return false;
    }

    const size_t spectrumSize = spectrum.size();
    // LOG_DEBUG("Initializing individual with spectrum size: {}", spectrumSize);  // Commented out - not essential

    // Create indices vector
    std::vector<int> indices(spectrumSize);
    std::iota(indices.begin(), indices.end(), 0);

    // Shuffle indices
    auto& rng = Random::instance();
    for (size_t i = spectrumSize; i > 1; --i) {
        size_t j = rng.getRandomSizeT(0, i - 1);
        std::swap(indices[i - 1], indices[j]);
    }

    // Set genes and validate
    individual.setGenes(indices);
    if (!validateGenes(indices, instance)) {
        // LOG_ERROR("Failed to validate initialized individual");  // Commented out - not essential
        return false;
    }

    LOG_DEBUG("Successfully initialized individual with {} genes", indices.size());
    return true;
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    
    if (!individual) {
        LOG_DEBUG("Validation failed: Null individual");
        return false;
    }
    
    // Check cache first
    auto it = m_validationCache.find(individual.get());
    if (it != m_validationCache.end()) {
        return it->second;
    }
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) {
        LOG_DEBUG("Validation failed: Empty genes vector in individual");
        m_validationCache[individual.get()] = false;
        return false;
    }
    
    bool result = validateGenes(genes, instance);
    m_validationCache[individual.get()] = result;
    
    // Clear cache if it gets too large
    if (m_validationCache.size() > 1000) {
        clearValidationCache();
    }
    
    return result;
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

std::string PermutationRepresentation::toDNA(const std::vector<size_t>& genes, const DNAInstance& instance) const {
    if (genes.empty()) return "";
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) return "";
    
    size_t k = instance.getK();
    size_t minLength = std::max(k, spectrum.size() / 6);  // More lenient minimum
    size_t maxLength = std::min(
        spectrum.size() * k / 2,     // Original calculation
        spectrum.size() + k * 2      // More reasonable upper bound
    );
    
    std::string assembled;
    assembled.reserve(maxLength);  // Pre-allocate to avoid reallocations
    
    size_t consecutivePoorOverlaps = 0;
    const size_t MAX_POOR_OVERLAPS = 8;
    
    // Start with first k-mer
    assembled = spectrum[genes[0]];
    
    for (size_t i = 1; i < genes.size(); ++i) {
        if (assembled.length() >= maxLength) {
            LOG_DEBUG("Stopping assembly - Length exceeded max: {}", assembled.length());
            break;  // Stop if we've exceeded max length
        }
        
        const std::string& current = spectrum[genes[i]];
        size_t bestOverlap = 0;
        size_t bestPosition = 0;
        bool foundGoodOverlap = false;
        
        // Try different overlap lengths
        for (size_t overlapLen = k-1; overlapLen >= k/4 && !foundGoodOverlap; --overlapLen) {
            // Check if adding this k-mer would exceed maxLength
            if (assembled.length() + current.length() - overlapLen > maxLength) {
                continue;  // Skip this overlap length if it would make sequence too long
            }
            
            size_t mismatches = countMismatches(assembled, current, overlapLen);
            if (mismatches <= k/3) {  // Allow more mismatches
                bestOverlap = overlapLen;
                bestPosition = assembled.length() - overlapLen;
                foundGoodOverlap = true;
                consecutivePoorOverlaps = 0;
                break;
            }
        }
        
        if (!foundGoodOverlap) {
            consecutivePoorOverlaps++;
            if (consecutivePoorOverlaps >= MAX_POOR_OVERLAPS) {
                LOG_DEBUG("Stopping assembly - Too many consecutive poor overlaps");
                break;
            }
            // Try minimal overlap if we can't find a good one
            if (assembled.length() + k - 1 <= maxLength) {
                assembled += current.substr(k-1);
            }
        } else {
            assembled += current.substr(bestOverlap);
        }
    }
    
    // Final length check
    if (assembled.length() < minLength || assembled.length() > maxLength) {
        LOG_DEBUG("Assembled sequence length invalid: {} (min: {}, max: {})",
                 assembled.length(), minLength, maxLength);
        return "";
    }
    
    return assembled;
}

bool PermutationRepresentation::validateGenes(
    const std::vector<size_t>& genes,
    const DNAInstance& instance) const {
    
    if (genes.empty()) return false;
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) return false;
    
    // Strict size check - must match spectrum size exactly for permutation
    if (genes.size() != spectrum.size()) {
        return false;
    }
    
    // Check for invalid indices and duplicates
    std::vector<bool> used(spectrum.size(), false);
    for (const auto& gene : genes) {
        if (gene >= spectrum.size()) return false;
        if (!instance.isRepAllowed() && used[gene]) return false;
        used[gene] = true;
    }
    
    // Generate DNA sequence
    std::string dna = toDNA(genes, instance);
    if (dna.empty()) return false;  // Length or coverage validation failed
    
    // Verify k-mer coverage
    size_t k = instance.getK();
    std::unordered_set<std::string> foundKmers;
    for (size_t i = 0; i <= dna.length() - k; ++i) {
        foundKmers.insert(dna.substr(i, k));
    }
    
    // More lenient coverage requirement - 8% minimum
    double coverage = static_cast<double>(foundKmers.size()) / spectrum.size();
    if (coverage < 0.08) return false;
    
    return true;
} 