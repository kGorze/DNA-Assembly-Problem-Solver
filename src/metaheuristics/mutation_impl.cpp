#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/utils/logging.h"
#include <memory>
#include <random>
#include <stdexcept>

PointMutation::PointMutation(double mutationRate) {
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        throw std::invalid_argument("Mutation rate must be between 0 and 1");
    }
    m_mutationRate = mutationRate;
    LOG_INFO("PointMutation initialized with rate: " + std::to_string(mutationRate));
}

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual) {
        LOG_ERROR("Null individual provided to mutation operator");
        return;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to mutation operator");
        return;
    }

    try {
        auto genes = individual->getGenes();
        if (genes.empty()) {
            LOG_ERROR("Empty genes vector in mutation");
            return;
        }
        if (genes.size() < 2) {
            LOG_WARNING("Individual has less than 2 genes, skipping mutation");
            return;
        }

        const auto& spectrum = instance.getSpectrum();
        if (spectrum.empty()) {
            LOG_ERROR("Empty spectrum in mutation");
            return;
        }

        const int realSpectrumSize = static_cast<int>(spectrum.size());
        const int k = instance.getK();
        if (k <= 1) {
            LOG_ERROR("Invalid k value in mutation: " + std::to_string(k));
            return;
        }

        // Validate all genes are within bounds and k-mers have sufficient length
        for (size_t i = 0; i < genes.size(); ++i) {
            if (genes[i] < 0 || genes[i] >= realSpectrumSize) {
                LOG_ERROR("Invalid gene index before mutation: " + std::to_string(genes[i]));
                return;
            }
            const auto& kmer = spectrum[genes[i]];
            if (kmer.length() < static_cast<size_t>(k)) {
                LOG_ERROR("K-mer too short at index " + std::to_string(i) + ": length=" + 
                         std::to_string(kmer.length()) + ", k=" + std::to_string(k));
                return;
            }
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        const int maxAttempts = 10;
        
        // Try mutation multiple times if needed
        for (int attempt = 0; attempt < maxAttempts; ++attempt) {
            auto mutatedGenes = genes;  // Reset to original genes for each attempt
            bool attemptMutated = false;
            
            // Attempt to mutate each gene with probability m_mutationRate
            for (size_t i = 0; i < genes.size(); ++i) {
                if (dis(gen) < m_mutationRate) {
                    // Find valid swap candidates that maintain k-mer connectivity
                    std::vector<size_t> validSwapPositions;
                    
                    // Check all possible swap positions
                    for (size_t j = 0; j < genes.size(); ++j) {
                        if (j == i) continue;
                        
                        // Try the swap
                        auto tempGenes = mutatedGenes;
                        std::swap(tempGenes[i], tempGenes[j]);
                        
                        // Check if the swap maintains valid k-mer sequence
                        bool isValid = true;
                        
                        // Check k-mers around position i
                        for (int offset = -1; offset <= 1 && isValid; ++offset) {
                            size_t pos = i + offset;
                            if (pos >= tempGenes.size() - 1) continue;
                            
                            const int curr_idx = tempGenes[pos];
                            const int next_idx = tempGenes[pos + 1];
                            
                            if (curr_idx < 0 || curr_idx >= realSpectrumSize ||
                                next_idx < 0 || next_idx >= realSpectrumSize) {
                                LOG_DEBUG("Invalid gene indices at positions " + std::to_string(pos) + 
                                        " and " + std::to_string(pos + 1));
                                isValid = false;
                                continue;
                            }
                            
                            const std::string& curr = spectrum[curr_idx];
                            const std::string& next = spectrum[next_idx];
                            
                            if (curr.length() < static_cast<size_t>(k) || next.length() < static_cast<size_t>(k)) {
                                LOG_DEBUG("K-mer too short at position " + std::to_string(pos) + 
                                        " or " + std::to_string(pos + 1));
                                isValid = false;
                                continue;
                            }
                            
                            try {
                                // Validate substring lengths before operation
                                if (curr.length() < static_cast<size_t>(k - 1)) {
                                    LOG_DEBUG("Current k-mer too short for suffix: length=" + 
                                            std::to_string(curr.length()) + ", need=" + 
                                            std::to_string(k - 1));
                                    isValid = false;
                                    continue;
                                }
                                if (next.length() < static_cast<size_t>(k - 1)) {
                                    LOG_DEBUG("Next k-mer too short for prefix: length=" + 
                                            std::to_string(next.length()) + ", need=" + 
                                            std::to_string(k - 1));
                                    isValid = false;
                                    continue;
                                }
                                
                                std::string suffix = curr.substr(curr.length() - (k-1));
                                std::string prefix = next.substr(0, k-1);
                                
                                int mismatches = 0;
                                for (int p = 0; p < k-1; ++p) {
                                    if (suffix[p] != prefix[p]) mismatches++;
                                }
                                if (mismatches > instance.getDeltaK()) {
                                    isValid = false;
                                }
                            } catch (const std::exception& e) {
                                LOG_ERROR("Substring error in mutation at position " + 
                                         std::to_string(pos) + ": " + std::string(e.what()));
                                isValid = false;
                                continue;
                            }
                        }
                        
                        // Check k-mers around position j (similar checks as above)
                        for (int offset = -1; offset <= 1 && isValid; ++offset) {
                            size_t pos = j + offset;
                            if (pos >= tempGenes.size() - 1) continue;
                            
                            const int curr_idx = tempGenes[pos];
                            const int next_idx = tempGenes[pos + 1];
                            
                            if (curr_idx < 0 || curr_idx >= realSpectrumSize ||
                                next_idx < 0 || next_idx >= realSpectrumSize) {
                                LOG_DEBUG("Invalid gene indices at positions " + std::to_string(pos) + 
                                        " and " + std::to_string(pos + 1));
                                isValid = false;
                                continue;
                            }
                            
                            const std::string& curr = spectrum[curr_idx];
                            const std::string& next = spectrum[next_idx];
                            
                            if (curr.length() < static_cast<size_t>(k) || next.length() < static_cast<size_t>(k)) {
                                LOG_DEBUG("K-mer too short at position " + std::to_string(pos) + 
                                        " or " + std::to_string(pos + 1));
                                isValid = false;
                                continue;
                            }
                            
                            try {
                                // Validate substring lengths before operation
                                if (curr.length() < static_cast<size_t>(k - 1)) {
                                    LOG_DEBUG("Current k-mer too short for suffix: length=" + 
                                            std::to_string(curr.length()) + ", need=" + 
                                            std::to_string(k - 1));
                                    isValid = false;
                                    continue;
                                }
                                if (next.length() < static_cast<size_t>(k - 1)) {
                                    LOG_DEBUG("Next k-mer too short for prefix: length=" + 
                                            std::to_string(next.length()) + ", need=" + 
                                            std::to_string(k - 1));
                                    isValid = false;
                                    continue;
                                }
                                
                                std::string suffix = curr.substr(curr.length() - (k-1));
                                std::string prefix = next.substr(0, k-1);
                                
                                int mismatches = 0;
                                for (int p = 0; p < k-1; ++p) {
                                    if (suffix[p] != prefix[p]) mismatches++;
                                }
                                if (mismatches > instance.getDeltaK()) {
                                    isValid = false;
                                }
                            } catch (const std::exception& e) {
                                LOG_ERROR("Substring error in mutation at position " + 
                                         std::to_string(pos) + ": " + std::string(e.what()));
                                isValid = false;
                                continue;
                            }
                        }
                        
                        if (isValid) {
                            validSwapPositions.push_back(j);
                        }
                    }
                    
                    // Handle empty validSwapPositions explicitly
                    if (validSwapPositions.empty()) {
                        LOG_DEBUG("No valid swap positions found for gene at position " + 
                                std::to_string(i));
                        continue;
                    }
                    
                    // Perform swap with a randomly chosen valid position
                    std::uniform_int_distribution<size_t> posDis(0, validSwapPositions.size() - 1);
                    size_t swapPos = validSwapPositions[posDis(gen)];
                    std::swap(mutatedGenes[i], mutatedGenes[swapPos]);
                    attemptMutated = true;
                }
            }
            
            // Only update the individual if the mutation was successful and produces a valid result
            if (attemptMutated) {
                // Validate all mutated genes are within bounds
                bool allValid = true;
                for (const auto& gene : mutatedGenes) {
                    if (gene < 0 || gene >= realSpectrumSize) {
                        LOG_ERROR("Invalid gene index after mutation: " + std::to_string(gene));
                        allValid = false;
                        break;
                    }
                }
                
                if (allValid) {
                    // Create new individual and validate it
                    auto mutatedIndividual = std::make_shared<Individual>(mutatedGenes);
                    if (representation->isValid(mutatedIndividual, instance)) {
                        // Double-check the mutated individual's genes
                        const auto& finalGenes = mutatedIndividual->getGenes();
                        if (!finalGenes.empty() && finalGenes.size() == genes.size()) {
                            individual = mutatedIndividual;
                            LOG_DEBUG("Mutation applied successfully");
                            return;  // Exit early on success
                        } else {
                            LOG_ERROR("Mutated individual has invalid gene count: " + 
                                     std::to_string(finalGenes.size()) + " vs original " + 
                                     std::to_string(genes.size()));
                        }
                    }
                }
            }
        }
        
        // If we get here, all mutation attempts failed
        LOG_WARNING("Failed to find valid mutation after " + std::to_string(maxAttempts) + 
                   " attempts - keeping original individual unchanged");
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error during mutation: " + std::string(e.what()));
    }
} 