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

class GuidedMutation : public IMutation {
private:
    double m_mutationRate;
    int m_minMutations;
    int m_maxAttempts;

public:
    explicit GuidedMutation(double mutationRate = 0.1, int minMutations = 1, int maxAttempts = 5)
        : m_mutationRate(mutationRate), m_minMutations(minMutations), m_maxAttempts(maxAttempts) {
        if (mutationRate < 0.0 || mutationRate > 1.0) {
            throw std::invalid_argument("Mutation rate must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

protected:
    // Helper methods for guided mutation
    bool tryReverseSegment(std::vector<int>& genes, size_t start, size_t length,
                          const DNAInstance& instance, std::shared_ptr<IRepresentation> representation);
    bool tryRealignSegment(std::vector<int>& genes, size_t start, size_t length,
                          const DNAInstance& instance, std::shared_ptr<IRepresentation> representation);
    bool tryMergeSubpaths(std::vector<int>& genes, size_t pos1, size_t pos2,
                         const DNAInstance& instance, std::shared_ptr<IRepresentation> representation);
};

class CombinedMutation : public IMutation {
private:
    std::shared_ptr<PointMutation> m_pointMutation;
    std::shared_ptr<SwapMutation> m_swapMutation;
    std::shared_ptr<GuidedMutation> m_guidedMutation;
    double m_swapProbability;
    double m_guidedProbability;
    double m_adaptiveMutationRate;
    int m_stagnationCounter;
    double m_lastBestFitness;

public:
    CombinedMutation(double pointRate = 0.1, double swapRate = 0.1, 
                     double swapProbability = 0.3, double guidedProbability = 0.4)
        : m_pointMutation(std::make_shared<PointMutation>(pointRate, 1))
        , m_swapMutation(std::make_shared<SwapMutation>(swapRate, 1))
        , m_guidedMutation(std::make_shared<GuidedMutation>(0.1, 1))
        , m_swapProbability(swapProbability)
        , m_guidedProbability(guidedProbability)
        , m_adaptiveMutationRate(pointRate)
        , m_stagnationCounter(0)
        , m_lastBestFitness(0.0) {
        if (swapProbability < 0.0 || swapProbability > 1.0) {
            throw std::invalid_argument("Swap probability must be between 0 and 1");
        }
        if (guidedProbability < 0.0 || guidedProbability > 1.0) {
            throw std::invalid_argument("Guided probability must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

    // Update mutation rate based on fitness improvement
    void updateMutationRate(double currentBestFitness);

private:
    static constexpr double MIN_MUTATION_RATE = 0.05;
    static constexpr double MAX_MUTATION_RATE = 0.2;
    static constexpr int STAGNATION_THRESHOLD = 5;
    static constexpr double IMPROVEMENT_THRESHOLD = 0.01;
}; 