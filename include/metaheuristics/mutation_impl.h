#pragma once

#include "../interfaces/i_mutation.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include "utils/random.h"
#include <random>
#include <algorithm>
#include <mutex>
#include "individual.h"
#include <memory>

class PointMutation : public IMutation {
private:
    double m_mutationRate;
    int m_minMutations;  // Minimum number of mutations to perform

public:
    explicit PointMutation(double mutationRate, int minMutations = 2) 
        : m_mutationRate(mutationRate), m_minMutations(minMutations) {
        if (mutationRate < 0.0 || mutationRate > 1.0) {
            throw std::invalid_argument("Mutation rate must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

    double getMutationRate() const { return m_mutationRate; }
};

class SwapMutation : public IMutation {
private:
    double m_mutationRate;
    int m_minSwaps;  // Minimum number of swaps to perform

public:
    explicit SwapMutation(double mutationRate, int minSwaps = 2) 
        : m_mutationRate(mutationRate), m_minSwaps(minSwaps) {
        if (mutationRate < 0.0 || mutationRate > 1.0) {
            throw std::invalid_argument("Mutation rate must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (!individual || !representation) {
            LOG_WARNING("Null individual or representation provided to mutation operator");
            return;
        }

        auto genes = individual->getGenes();
        if (genes.size() < 2) {
            LOG_WARNING("Individual has less than 2 genes, cannot perform mutation");
            return;
        }

        // Create a copy for mutation
        auto mutated = std::make_shared<Individual>(genes);
        auto& mutatedGenes = mutated->getGenes();
        
        // Calculate number of swaps based on mutation rate and size, ensuring minimum swaps
        int numSwaps = std::max(m_minSwaps, static_cast<int>(genes.size() * m_mutationRate));
        
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
                int pos1 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
                int pos2;
                do {
                    pos2 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
                } while (pos1 == pos2);
                
                // Perform swap
                std::swap(mutatedGenes[pos1], mutatedGenes[pos2]);
                
                // Check if this swap made any improvement
                auto tempIndividual = std::make_shared<Individual>(mutatedGenes);
                if (representation->isValid(tempIndividual, instance)) {
                    swapSuccessful = true;
                    anyValidMutation = true;
                    consecutiveFailures = 0;
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
};

class CombinedMutation : public IMutation {
private:
    std::shared_ptr<PointMutation> m_pointMutation;
    std::shared_ptr<SwapMutation> m_swapMutation;
    double m_swapProbability;
    double m_bothProbability;

public:
    CombinedMutation(double pointRate, double swapRate, double swapProbability = 0.5, double bothProbability = 0.4) 
        : m_pointMutation(std::make_shared<PointMutation>(pointRate, 2))
        , m_swapMutation(std::make_shared<SwapMutation>(swapRate, 2))
        , m_swapProbability(swapProbability)
        , m_bothProbability(bothProbability) {
        if (swapProbability < 0.0 || swapProbability > 1.0) {
            throw std::invalid_argument("Swap probability must be between 0 and 1");
        }
        if (bothProbability < 0.0 || bothProbability > 1.0) {
            throw std::invalid_argument("Both probability must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
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
            
            // With bothProbability chance, apply both mutations
            if (rand < m_bothProbability) {
                m_swapMutation->mutate(currentMutated, instance, representation);
                m_pointMutation->mutate(currentMutated, instance, representation);
            }
            // Otherwise choose one based on swapProbability
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
                LOG_DEBUG("CombinedMutation: Successfully applied mutation on attempt " + std::to_string(attempt + 1));
            }
        }
        
        // Update individual if mutation was successful
        if (mutationSuccessful) {
            individual = mutated;
        } else {
            LOG_DEBUG("CombinedMutation: Failed to find valid mutation after 3 attempts");
        }
    }
}; 