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
        LOG_ERROR("Cannot perform point mutation: null individual or representation");
        return;
    }

    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("Cannot perform point mutation: spectrum is empty");
        return;
    }

    const size_t spectrumSize = spectrum.size();
    auto genes = individual->getGenes();

    if (genes.empty()) {
        LOG_ERROR("Cannot perform point mutation: empty gene vector");
        return;
    }

    // Validate all genes are within spectrum size
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] < 0 || static_cast<size_t>(genes[i]) >= spectrumSize) {
            LOG_ERROR("Invalid gene value at position {}: {} (spectrum size: {})", 
                     i, genes[i], spectrumSize);
            return;
        }
    }

    auto& rng = Random::instance();
    const int k = instance.getK(); // Get k-mer size for overlap validation
    
    // Calculate maximum mutations based on sequence length and mutation rate
    const size_t maxPossibleMutations = std::min(
        static_cast<size_t>(genes.size() * m_mutationRate),
        static_cast<size_t>(3)  // Hard cap at 3 mutations
    );
    
    if (maxPossibleMutations == 0) {
        LOG_ERROR("Cannot determine valid number of mutations");
        return;
    }

    const size_t numMutations = rng.getRandomSizeT(1, maxPossibleMutations);
    LOG_DEBUG("Attempting {} mutations on individual with {} genes", numMutations, genes.size());
    
    size_t successfulMutations = 0;
    for (size_t i = 0; i < numMutations; ++i) {
        // Try multiple positions to find a good mutation
        const size_t maxAttempts = 10;
        bool foundValidMutation = false;
        
        for (size_t attempt = 0; attempt < maxAttempts; ++attempt) {
            // Generate position with strict bounds checking
            const size_t pos = rng.getRandomSizeT(0, genes.size() - 1);
            if (pos >= genes.size()) {
                continue;
            }

            const int originalValue = genes[pos];
            std::vector<int> candidates;
            candidates.reserve(spectrumSize);

            // Build list of candidate values with valid overlaps
            for (size_t newValue = 0; newValue < spectrumSize; ++newValue) {
                if (static_cast<int>(newValue) == originalValue) {
                    continue;
                }

                genes[pos] = static_cast<int>(newValue);
                bool validOverlaps = true;

                // Check overlaps with neighbors
                if (pos > 0) {
                    double overlapBefore = dna_utils::calculateEdgeWeight(
                        spectrum[genes[pos-1]], spectrum[genes[pos]], k);
                    if (overlapBefore <= 0) {
                        validOverlaps = false;
                    }
                }
                if (pos < genes.size() - 1 && validOverlaps) {
                    double overlapAfter = dna_utils::calculateEdgeWeight(
                        spectrum[genes[pos]], spectrum[genes[pos+1]], k);
                    if (overlapAfter <= 0) {
                        validOverlaps = false;
                    }
                }

                if (validOverlaps) {
                    candidates.push_back(static_cast<int>(newValue));
                }
            }

            // Restore original value
            genes[pos] = originalValue;

            if (!candidates.empty()) {
                // Choose a random candidate
                const size_t candidateIdx = rng.getRandomSizeT(0, candidates.size() - 1);
                genes[pos] = candidates[candidateIdx];

                // Validate the complete sequence
                auto tempIndividual = std::make_shared<Individual>(genes);
                if (representation->isValid(tempIndividual, instance)) {
                    foundValidMutation = true;
                    successfulMutations++;
                    break;
                }

                // Revert if not valid
                genes[pos] = originalValue;
            }
        }

        if (!foundValidMutation) {
            LOG_DEBUG("Failed to find valid mutation after {} attempts", maxAttempts);
        }
    }

    if (successfulMutations > 0) {
        individual->setGenes(genes);
        LOG_DEBUG("Successfully performed {} out of {} attempted mutations", 
                 successfulMutations, numMutations);
    } else {
        LOG_WARNING("No successful mutations performed");
    }
}

void SwapMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual || !representation) {
        LOG_ERROR("Cannot perform swap mutation: null individual or representation");
        return;
    }

    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("Cannot perform swap mutation: spectrum is empty");
        return;
    }

    const size_t spectrumSize = spectrum.size();
    auto genes = individual->getGenes();

    // Strict validation of gene vector size
    if (genes.empty() || genes.size() < 2) {
        LOG_ERROR("Cannot perform swap mutation: need at least 2 genes");
        return;
    }

    if (genes.size() != spectrumSize) {
        LOG_ERROR("Gene vector size ({}) does not match spectrum size ({})", genes.size(), spectrumSize);
        return;
    }

    // Validate all genes are within spectrum size
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] < 0 || static_cast<size_t>(genes[i]) >= spectrumSize) {
            LOG_ERROR("Invalid gene value at position {}: {} (spectrum size: {}). Skipping mutation.", 
                     i, genes[i], spectrumSize);
            return;
        }
    }

    auto& rng = Random::instance();
    // Limit number of swaps to a reasonable value based on sequence length
    const size_t maxSwaps = std::min(size_t{3}, size_t{genes.size() / 8});
    const size_t numSwaps = rng.getRandomSizeT(1, maxSwaps);

    LOG_DEBUG("Attempting {} swaps on individual with {} genes", numSwaps, genes.size());

    size_t successfulSwaps = 0;
    const int k = instance.getK(); // Get k-mer size for overlap validation
    
    for (size_t i = 0; i < numSwaps; ++i) {
        // Try multiple positions to find a good swap
        const size_t maxAttempts = 10;
        bool foundValidSwap = false;
        
        for (size_t attempt = 0; attempt < maxAttempts; ++attempt) {
            // Generate positions with strict bounds checking
            const size_t pos1 = rng.getRandomSizeT(0, genes.size() - 1);
            const size_t pos2 = rng.getRandomSizeT(0, genes.size() - 1);
            
            if (pos1 == pos2 || pos1 >= genes.size() || pos2 >= genes.size()) {
                continue;
            }

            // Store original values
            [[maybe_unused]] const int gene1 = genes[pos1];
            [[maybe_unused]] const int gene2 = genes[pos2];

            // Perform the swap
            std::swap(genes[pos1], genes[pos2]);

            // Check k-mer overlaps around the swapped positions
            bool validOverlaps = true;
            for (size_t pos : {pos1, pos2}) {
                if (pos > 0) {
                    double overlapBefore = dna_utils::calculateEdgeWeight(
                        spectrum[genes[pos-1]], spectrum[genes[pos]], k);
                    if (overlapBefore <= 0) {
                        validOverlaps = false;
                        break;
                    }
                }
                if (pos < genes.size() - 1) {
                    double overlapAfter = dna_utils::calculateEdgeWeight(
                        spectrum[genes[pos]], spectrum[genes[pos+1]], k);
                    if (overlapAfter <= 0) {
                        validOverlaps = false;
                        break;
                    }
                }
            }

            // If overlaps are invalid, try to repair
            if (!validOverlaps) {
                // Try to find a better position for gene2 near pos1
                bool repaired = false;
                for (int offset = -2; offset <= 2 && !repaired; ++offset) {
                    size_t newPos = pos1;
                    if (offset < 0 && static_cast<int>(pos1) + offset >= 0) {
                        newPos = pos1 + offset;
                    } else if (offset > 0 && pos1 + offset < genes.size()) {
                        newPos = pos1 + offset;
                    } else {
                        continue;
                    }
                    
                    // Try placing gene2 at the new position
                    std::swap(genes[newPos], genes[pos2]);
                    
                    // Validate overlaps at both positions
                    bool validRepair = true;
                    for (size_t p : {newPos, pos2}) {
                        if (p > 0) {
                            double overlap = dna_utils::calculateEdgeWeight(
                                spectrum[genes[p-1]], spectrum[genes[p]], k);
                            if (overlap <= 0) {
                                validRepair = false;
                                break;
                            }
                        }
                        if (p < genes.size() - 1) {
                            double overlap = dna_utils::calculateEdgeWeight(
                                spectrum[genes[p]], spectrum[genes[p+1]], k);
                            if (overlap <= 0) {
                                validRepair = false;
                                break;
                            }
                        }
                    }
                    
                    if (validRepair) {
                        repaired = true;
                        validOverlaps = true;
                    } else {
                        // Undo the repair attempt
                        std::swap(genes[newPos], genes[pos2]);
                    }
                }
            }

            if (validOverlaps) {
                // Create temporary individual to validate the complete sequence
                auto tempIndividual = std::make_shared<Individual>(genes);
                if (representation->isValid(tempIndividual, instance)) {
                    foundValidSwap = true;
                    successfulSwaps++;
                    break;
                }
            }

            // Revert the swap if not valid
            std::swap(genes[pos1], genes[pos2]);
        }

        if (!foundValidSwap) {
            LOG_DEBUG("Failed to find valid swap after {} attempts", maxAttempts);
        }
    }

    if (successfulSwaps > 0) {
        individual->setGenes(genes);
        LOG_DEBUG("Successfully performed {} out of {} attempted swaps", 
                 successfulSwaps, numSwaps);
    } else {
        LOG_WARNING("No successful swaps performed during mutation");
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