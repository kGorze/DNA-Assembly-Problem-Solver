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
        if (genes.size() < 2) {
            LOG_WARNING("Individual has less than 2 genes, skipping mutation");
            return;
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        bool mutated = false;
        std::vector<int> mutatedGenes = genes;
        
        // Calculate real spectrum size
        const int realSpectrumSize = instance.getN() - instance.getK() + 1;
        
        // Higher chance to mutate error k-mers
        for (size_t i = 0; i < genes.size(); ++i) {
            double mutationChance = m_mutationRate;
            
            // Increase mutation chance for error k-mers
            if (genes[i] >= realSpectrumSize) {
                mutationChance *= 2.0;  // Double the mutation rate for error k-mers
            }
            
            if (dis(gen) < mutationChance) {
                // Prefer swapping with real k-mers when possible
                std::vector<size_t> validSwapPositions;
                for (size_t j = 0; j < genes.size(); ++j) {
                    if (j != i) {
                        // Prioritize real k-mers as swap targets
                        if (genes[j] < realSpectrumSize) {
                            validSwapPositions.push_back(j);
                            validSwapPositions.push_back(j);  // Add twice to increase probability
                        } else {
                            validSwapPositions.push_back(j);
                        }
                    }
                }
                
                if (!validSwapPositions.empty()) {
                    std::uniform_int_distribution<size_t> posDis(0, validSwapPositions.size() - 1);
                    size_t swapPos = validSwapPositions[posDis(gen)];
                    std::swap(mutatedGenes[i], mutatedGenes[swapPos]);
                    mutated = true;
                }
            }
        }

        if (mutated) {
            auto mutatedIndividual = std::make_shared<Individual>(mutatedGenes);
            if (representation->isValid(mutatedIndividual, instance)) {
                individual = mutatedIndividual;
                LOG_DEBUG("Mutation applied successfully");
            } else {
                LOG_WARNING("Mutation produced invalid individual, keeping original");
            }
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Error during mutation: " + std::string(e.what()));
    }
} 