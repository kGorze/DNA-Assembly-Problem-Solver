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
    
    if (!individual) {
        LOG_ERROR("Individual is null");
        return;
    }

    auto genes = individual->getGenes();
    if (genes.size() < 2) {
        LOG_ERROR("Individual has less than 2 genes, cannot perform mutation");
        return;
    }

    // Create a copy for mutation
    auto mutated = std::make_shared<Individual>(genes);
    auto& mutatedGenes = mutated->getGenes();
    
    // Calculate number of mutations based on size, ensuring minimum mutations
    int numMutations = std::max(m_minMutations, static_cast<int>(genes.size() * m_mutationRate));
    
    auto& rng = Random::instance();
    bool anyValidMutation = false;
    int consecutiveFailures = 0;
    
    // Try multiple mutations with backtracking
    for (int i = 0; i < numMutations && consecutiveFailures < 5; ++i) {
        // Store current state
        auto currentState = mutatedGenes;
        
        // Try up to 3 different positions for a successful mutation
        bool mutationSuccessful = false;
        for (int attempt = 0; attempt < 3 && !mutationSuccessful; ++attempt) {
            // Select two random positions for swapping
            int pos1 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
            int pos2;
            do {
                pos2 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
            } while (pos1 == pos2);
            
            // Perform swap
            std::swap(mutatedGenes[pos1], mutatedGenes[pos2]);
            
            // Try additional random swap to increase diversity
            if (rng.generateProbability() < 0.3) {  // 30% chance for additional swap
                int pos3 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
                int pos4;
                do {
                    pos4 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
                } while (pos3 == pos4);
                std::swap(mutatedGenes[pos3], mutatedGenes[pos4]);
            }
            
            // Check if this mutation made any improvement
            auto tempIndividual = std::make_shared<Individual>(mutatedGenes);
            if (representation->isValid(tempIndividual, instance)) {
                mutationSuccessful = true;
                anyValidMutation = true;
                consecutiveFailures = 0;
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
    
    // Only update if we made valid changes and the result is different
    if (anyValidMutation && mutatedGenes != genes) {
        individual = mutated;
        LOG_DEBUG("PointMutation: Successfully performed multiple mutations");
    } else {
        LOG_DEBUG("PointMutation: No valid mutations found or no changes made");
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

    auto genes = individual->getGenes();
    if (genes.size() < 4) {  // Need at least 4 genes for meaningful operations
        return;
    }

    auto& rng = Random::instance();
    bool mutationSuccessful = false;
    std::vector<int> bestGenes = genes;
    double bestQuality = -1.0;

    // Try different mutation strategies
    for (int attempt = 0; attempt < m_maxAttempts && !mutationSuccessful; ++attempt) {
        auto currentGenes = genes;
        bool improved = false;

        // Choose mutation strategy based on probabilities
        double strategy = rng.generateProbability();
        
        if (strategy < 0.4) {  // 40% chance for segment reversal
            size_t start = rng.getRandomInt(0, static_cast<int>(genes.size() - 4));
            size_t length = rng.getRandomInt(3, static_cast<int>(std::min(genes.size() - start, size_t(8))));
            improved = tryReverseSegment(currentGenes, start, length, instance, representation);
        }
        else if (strategy < 0.7) {  // 30% chance for segment realignment
            size_t start = rng.getRandomInt(0, static_cast<int>(genes.size() - 4));
            size_t length = rng.getRandomInt(3, static_cast<int>(std::min(genes.size() - start, size_t(6))));
            improved = tryRealignSegment(currentGenes, start, length, instance, representation);
        }
        else {  // 30% chance for subpath merging
            size_t pos1 = rng.getRandomInt(0, static_cast<int>(genes.size() - 3));
            size_t pos2 = rng.getRandomInt(0, static_cast<int>(genes.size() - 3));
            if (pos1 > pos2) std::swap(pos1, pos2);
            improved = tryMergeSubpaths(currentGenes, pos1, pos2, instance, representation);
        }

        if (improved) {
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
                }
            }
        }
    }

    if (mutationSuccessful) {
        individual = std::make_shared<Individual>(bestGenes);
        LOG_DEBUG("GuidedMutation: Successfully applied mutation");
    }
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
    
    // Try to merge subpaths by finding best connection point
    size_t bestPos = pos1;
    double bestQuality = -1.0;
    
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
                bestPos = i;
                genes = merged;
            }
        }
    }
    
    if (bestQuality > -1.0) {
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
        LOG_WARNING("Null individual or representation in combined mutation");
        return;
    }

    // Create a copy for mutation
    auto mutated = std::make_shared<Individual>(individual->getGenes());
    bool mutationSuccessful = false;
    
    // Try up to 3 times to get a valid mutation
    for (int attempt = 0; attempt < 3 && !mutationSuccessful; attempt++) {
        auto currentMutated = std::make_shared<Individual>(individual->getGenes());
        
        auto& rng = Random::instance();
        double rand = rng.generateProbability();
        
        // Apply guided mutation with guidedProbability chance
        if (rand < m_guidedProbability) {
            m_guidedMutation->mutate(currentMutated, instance, representation);
        }
        // Otherwise choose between point and swap mutation
        else {
            if (rng.generateProbability() < m_swapProbability) {
                m_swapMutation->mutate(currentMutated, instance, representation);
            } else {
                m_pointMutation->mutate(currentMutated, instance, representation);
            }
        }
        
        // Check if mutation was successful and different from original
        if (representation->isValid(currentMutated, instance) && 
            currentMutated->getGenes() != individual->getGenes()) {
            mutated = currentMutated;
            mutationSuccessful = true;
            LOG_DEBUG("CombinedMutation: Successfully applied mutation on attempt " + 
                     std::to_string(attempt + 1));
        }
    }
    
    // Update individual if mutation was successful
    if (mutationSuccessful) {
        individual = mutated;
    } else {
        LOG_DEBUG("CombinedMutation: Failed to find valid mutation after 3 attempts");
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
            m_stagnationCounter = 0;
        }
    }
    
    // Update mutation rates for component operators
    m_pointMutation = std::make_shared<PointMutation>(m_adaptiveMutationRate, 1);
    m_swapMutation = std::make_shared<SwapMutation>(m_adaptiveMutationRate, 1);
    m_guidedMutation = std::make_shared<GuidedMutation>(m_adaptiveMutationRate, 1);
    
    m_lastBestFitness = currentBestFitness;
    
    LOG_DEBUG("Updated mutation rate to: " + std::to_string(m_adaptiveMutationRate));
} 