#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/dna_utils.h"
#include <memory>
#include <random>
#include <stdexcept>
#include <sstream>
#include <algorithm>

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual || !representation) {
        LOG_WARNING("Null individual or representation provided to mutation operator");
        return;
    }

    auto genes = individual->getGenes();
    const size_t spectrumSize = instance.getSpectrum().size();
    
    // Validate gene vector size
    if (genes.size() != spectrumSize) {
        LOG_WARNING("Individual's genes size (" + std::to_string(genes.size()) + 
                   ") does not match spectrum size (" + std::to_string(spectrumSize) + ")");
        return;
    }

    // Create a copy for mutation
    auto mutated = std::make_shared<Individual>(genes);
    auto& mutatedGenes = mutated->getGenes();
    
    // Calculate number of mutations based on mutation rate
    int numMutations = std::max(1, static_cast<int>(spectrumSize * m_mutationRate));
    
    auto& rng = Random::instance();
    bool anyValidMutation = false;
    int consecutiveFailures = 0;
    
    // Try multiple point mutations with backtracking
    for (int i = 0; i < numMutations && consecutiveFailures < 5; ++i) {
        // Store current state
        auto currentState = mutatedGenes;
        
        // Try up to 3 different positions for a successful mutation
        bool mutationSuccessful = false;
        for (int attempt = 0; attempt < 3 && !mutationSuccessful; ++attempt) {
            int pos = rng.getRandomInt(0, static_cast<int>(spectrumSize - 1));
            
            // Generate a new value that's different from the current one and within bounds
            int currentValue = mutatedGenes[pos];
            int newValue;
            do {
                newValue = rng.getRandomInt(0, static_cast<int>(spectrumSize - 1));
            } while (newValue == currentValue || newValue < 0 || newValue >= static_cast<int>(spectrumSize));
            
            // Validate indices before mutation
            if (pos < 0 || pos >= static_cast<int>(spectrumSize)) {
                LOG_WARNING("Invalid position index in point mutation");
                continue;
            }
            
            mutatedGenes[pos] = newValue;
            
            // Validate all genes after mutation
            bool validIndices = true;
            for (int gene : mutatedGenes) {
                if (gene < 0 || gene >= static_cast<int>(spectrumSize)) {
                    validIndices = false;
                    break;
                }
            }
            
            if (!validIndices) {
                mutatedGenes = currentState;
                continue;
            }
            
            // Check if this mutation made any improvement
            auto tempIndividual = std::make_shared<Individual>(mutatedGenes);
            if (representation->isValid(tempIndividual, instance)) {
                mutationSuccessful = true;
                anyValidMutation = true;
                consecutiveFailures = 0;
                LOG_DEBUG("PointMutation: Successfully mutated position " + std::to_string(pos) + 
                         " from " + std::to_string(currentValue) + " to " + std::to_string(newValue));
            } else {
                // Undo this mutation and try another position
                mutatedGenes = currentState;
            }
        }
        
        if (!mutationSuccessful) {
            consecutiveFailures++;
            mutatedGenes = currentState;  // Revert to last valid state
        }
    }
    
    // Only update if we made valid changes
    if (anyValidMutation) {
        individual = mutated;
        LOG_DEBUG("PointMutation: Successfully performed multiple mutations");
    }
}

void SwapMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual || !representation) {
        LOG_WARNING("Null individual or representation provided to mutation operator");
        return;
    }

    auto genes = individual->getGenes();
    const size_t spectrumSize = instance.getSpectrum().size();
    
    // Validate gene vector size
    if (genes.size() != spectrumSize) {
        LOG_WARNING("Individual's genes size (" + std::to_string(genes.size()) + 
                   ") does not match spectrum size (" + std::to_string(spectrumSize) + ")");
        return;
    }

    // Create a copy for mutation
    auto mutated = std::make_shared<Individual>(genes);
    auto& mutatedGenes = mutated->getGenes();
    
    // Calculate number of swaps based on mutation rate
    int numSwaps = std::max(1, static_cast<int>(spectrumSize * m_mutationRate));
    
    auto& rng = Random::instance();
    bool anyValidMutation = false;
    int consecutiveFailures = 0;
    
    // Try multiple swaps with backtracking
    for (int i = 0; i < numSwaps && consecutiveFailures < 5; ++i) {
        // Store current state
        auto currentState = mutatedGenes;
        
        // Try up to 3 different positions for a successful swap
        bool swapSuccessful = false;
        for (int attempt = 0; attempt < 3 && !swapSuccessful; ++attempt) {
            int pos1 = rng.getRandomInt(0, static_cast<int>(spectrumSize - 1));
            int pos2;
            do {
                pos2 = rng.getRandomInt(0, static_cast<int>(spectrumSize - 1));
            } while (pos1 == pos2);
            
            // Validate indices before swap
            if (pos1 < 0 || pos1 >= static_cast<int>(spectrumSize) || 
                pos2 < 0 || pos2 >= static_cast<int>(spectrumSize)) {
                LOG_WARNING("Invalid position indices in swap mutation");
                continue;
            }
            
            // Perform swap
            std::swap(mutatedGenes[pos1], mutatedGenes[pos2]);
            
            // Validate all genes after swap
            bool validIndices = true;
            for (int gene : mutatedGenes) {
                if (gene < 0 || gene >= static_cast<int>(spectrumSize)) {
                    validIndices = false;
                    break;
                }
            }
            
            if (!validIndices) {
                mutatedGenes = currentState;
                continue;
            }
            
            // Check if this swap made any improvement
            auto tempIndividual = std::make_shared<Individual>(mutatedGenes);
            if (representation->isValid(tempIndividual, instance)) {
                swapSuccessful = true;
                anyValidMutation = true;
                consecutiveFailures = 0;
                LOG_DEBUG("SwapMutation: Successfully swapped positions " + std::to_string(pos1) + 
                         " and " + std::to_string(pos2));
            } else {
                // Undo this swap and try another position
                mutatedGenes = currentState;
            }
        }
        
        if (!swapSuccessful) {
            consecutiveFailures++;
            mutatedGenes = currentState;  // Revert to last valid state
        }
    }
    
    // Only update if we made valid changes
    if (anyValidMutation) {
        individual = mutated;
        LOG_DEBUG("SwapMutation: Successfully performed multiple swaps");
    }
}

void GuidedMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual || !representation) {
        LOG_WARNING("Null individual or representation in guided mutation");
        return;
    }

    auto originalGenes = individual->getGenes();
    if (originalGenes.size() < 4) {  // Need at least 4 genes for meaningful operations
        return;
    }

    const size_t spectrumSize = instance.getSpectrum().size();
    
    // Validate original genes first
    bool hasInvalidGenes = false;
    for (int gene : originalGenes) {
        if (gene < 0 || gene >= static_cast<int>(spectrumSize)) {
            hasInvalidGenes = true;
            break;
        }
    }
    if (hasInvalidGenes) {
        LOG_WARNING("Original genes contain invalid values, skipping guided mutation");
        return;
    }

    auto& rng = Random::instance();
    bool mutationSuccessful = false;
    std::vector<int> bestGenes = originalGenes;
    double bestQuality = -1.0;

    // Try different mutation strategies
    for (int attempt = 0; attempt < m_maxAttempts && !mutationSuccessful; ++attempt) {
        auto currentGenes = originalGenes;
        bool improved = false;

        // Choose mutation strategy based on probabilities
        double strategy = rng.generateProbability();
        
        if (strategy < 0.4) {  // 40% chance for segment reversal
            size_t start = rng.getRandomInt(0, static_cast<int>(originalGenes.size() - 4));
            size_t length = rng.getRandomInt(3, static_cast<int>(std::min(originalGenes.size() - start, size_t(8))));
            improved = tryReverseSegment(currentGenes, start, length, instance, representation);
        }
        else if (strategy < 0.7) {  // 30% chance for segment realignment
            size_t start = rng.getRandomInt(0, static_cast<int>(originalGenes.size() - 4));
            size_t length = rng.getRandomInt(3, static_cast<int>(std::min(originalGenes.size() - start, size_t(6))));
            improved = tryRealignSegment(currentGenes, start, length, instance, representation);
        }
        else {  // 30% chance for subpath merging
            size_t pos1 = rng.getRandomInt(0, static_cast<int>(originalGenes.size() - 3));
            size_t pos2 = rng.getRandomInt(0, static_cast<int>(originalGenes.size() - 3));
            if (pos1 > pos2) std::swap(pos1, pos2);
            improved = tryMergeSubpaths(currentGenes, pos1, pos2, instance, representation);
        }

        if (improved) {
            // Validate all genes after mutation
            bool validIndices = true;
            for (int gene : currentGenes) {
                if (gene < 0 || gene >= static_cast<int>(spectrumSize)) {
                    validIndices = false;
                    break;
                }
            }

            if (validIndices) {
                auto tempIndividual = std::make_shared<Individual>(currentGenes);
                if (representation->isValid(tempIndividual, instance)) {
                    // Calculate quality of the mutation
                    double quality = 0.0;
                    const auto& spectrum = instance.getSpectrum();
                    for (size_t i = 0; i < currentGenes.size() - 1; ++i) {
                        quality += dna_utils::calculateEdgeWeight(
                            spectrum[currentGenes[i]], 
                            spectrum[currentGenes[i + 1]], 
                            instance.getK());
                    }

                    if (quality > bestQuality) {
                        bestQuality = quality;
                        bestGenes = currentGenes;
                        mutationSuccessful = true;
                        LOG_DEBUG("GuidedMutation: Found better solution with quality " + 
                                 std::to_string(quality));
                    }
                }
            }
        }
    }

    // If no successful mutation was found, keep the original genes
    if (!mutationSuccessful) {
        LOG_DEBUG("GuidedMutation: No successful mutation found, keeping original genes");
        return;
    }

    individual->setGenes(bestGenes);
    LOG_DEBUG("GuidedMutation: Successfully performed mutation");
}

bool GuidedMutation::tryReverseSegment(
    std::vector<int>& genes,
    size_t start,
    size_t length,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (start + length > genes.size()) return false;

    // Store original segment
    std::vector<int> originalSegment(genes.begin() + start, genes.begin() + start + length);
    
    // Reverse the segment
    std::reverse(genes.begin() + start, genes.begin() + start + length);
    
    // Check if reversal improves overlap quality
    auto tempIndividual = std::make_shared<Individual>(genes);
    if (!representation->isValid(tempIndividual, instance)) {
        // Restore original segment if invalid
        std::copy(originalSegment.begin(), originalSegment.end(), genes.begin() + start);
        return false;
    }
    
    return true;
}

bool GuidedMutation::tryRealignSegment(
    std::vector<int>& genes,
    size_t start,
    size_t length,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (start + length > genes.size()) return false;

    // Store original segment
    std::vector<int> originalSegment(genes.begin() + start, genes.begin() + start + length);
    
    // Try different alignments of the segment
    std::vector<int> bestAlignment = originalSegment;
    double bestQuality = -1.0;
    
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    // Try each possible rotation of the segment
    for (size_t i = 0; i < length; ++i) {
        std::rotate(genes.begin() + start, genes.begin() + start + 1, genes.begin() + start + length);
        
        // Calculate quality of this alignment
        double quality = 0.0;
        for (size_t j = start; j < start + length - 1; ++j) {
            quality += dna_utils::calculateEdgeWeight(
                spectrum[genes[j]], spectrum[genes[j + 1]], k);
        }
        
        if (quality > bestQuality) {
            auto tempIndividual = std::make_shared<Individual>(genes);
            if (representation->isValid(tempIndividual, instance)) {
                bestQuality = quality;
                bestAlignment = std::vector<int>(
                    genes.begin() + start, 
                    genes.begin() + start + length);
            }
        }
    }
    
    if (bestQuality > -1.0) {
        std::copy(bestAlignment.begin(), bestAlignment.end(), genes.begin() + start);
        return true;
    }
    
    // Restore original segment if no improvement found
    std::copy(originalSegment.begin(), originalSegment.end(), genes.begin() + start);
    return false;
}

bool GuidedMutation::tryMergeSubpaths(
    std::vector<int>& genes,
    size_t pos1,
    size_t pos2,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (pos1 >= pos2 || pos2 >= genes.size() - 2) return false;

    // Store original sequence
    std::vector<int> original = genes;
    
    double bestQuality = -1.0;
    std::vector<int> bestMerge = original;
    
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    // Try different merge points
    for (size_t i = pos1; i <= pos2; ++i) {
        // Temporarily merge at this position
        std::vector<int> merged;
        merged.insert(merged.end(), genes.begin(), genes.begin() + i);
        merged.insert(merged.end(), genes.begin() + pos2, genes.begin() + pos2 + 2);
        merged.insert(merged.end(), genes.begin() + i, genes.begin() + pos2);
        merged.insert(merged.end(), genes.begin() + pos2 + 2, genes.end());
        
        // Calculate quality of this merge
        double quality = 0.0;
        for (size_t j = 0; j < merged.size() - 1; ++j) {
            quality += dna_utils::calculateEdgeWeight(
                spectrum[merged[j]], spectrum[merged[j + 1]], k);
        }
        
        if (quality > bestQuality) {
            auto tempIndividual = std::make_shared<Individual>(merged);
            if (representation->isValid(tempIndividual, instance)) {
                bestQuality = quality;
                bestMerge = merged;
            }
        }
    }
    
    if (bestQuality > -1.0) {
        genes = bestMerge;
        return true;
    }
    
    // Restore original if no improvement found
    genes = original;
    return false;
}

void CombinedMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual || !representation) {
        LOG_WARNING("Null individual or representation provided to mutation operator");
        return;
    }

    auto& rng = Random::instance();
    double roll = rng.generateProbability();

    // 40% chance for guided mutation
    if (roll < 0.4 && m_guidedMutation) {
        m_guidedMutation->mutate(individual, instance, representation);
    }
    // 30% chance for swap mutation
    else if (roll < 0.7 && m_swapMutation) {
        m_swapMutation->mutate(individual, instance, representation);
    }
    // 30% chance for point mutation
    else if (m_pointMutation) {
        m_pointMutation->mutate(individual, instance, representation);
    }
}

void CombinedMutation::updateMutationRate(double currentBestFitness) {
    // Calculate improvement
    double improvement = currentBestFitness - m_lastBestFitness;
    
    if (improvement > IMPROVEMENT_THRESHOLD) {
        // Reset stagnation counter and slightly decrease mutation rate
        m_stagnationCounter = 0;
        m_adaptiveMutationRate = std::max(
            MIN_MUTATION_RATE,
            m_adaptiveMutationRate * 0.95
        );
    } else {
        // Increment stagnation counter
        m_stagnationCounter++;
        
        // If stagnated for too long, increase mutation rate
        if (m_stagnationCounter >= STAGNATION_THRESHOLD) {
            m_adaptiveMutationRate = std::min(
                MAX_MUTATION_RATE,
                m_adaptiveMutationRate * 1.1
            );
            
            // Reset counter to avoid too rapid increase
            m_stagnationCounter = 0;
        }
    }
    
    // Update mutation rate and propagate to component mutations
    setMutationRate(m_adaptiveMutationRate);
    m_lastBestFitness = currentBestFitness;
} 