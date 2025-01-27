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

public:
    explicit PointMutation(double mutationRate = 0.1) : m_mutationRate(mutationRate) {
        if (mutationRate < 0.0 || mutationRate > 1.0) {
            throw std::invalid_argument("Mutation rate must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

    void setMutationRate(double rate) override {
        m_mutationRate = std::clamp(rate, 0.1, 0.4);
    }

    double getMutationRate() const override {
        return m_mutationRate;
    }
};

class SwapMutation : public IMutation {
private:
    double m_mutationRate;

public:
    explicit SwapMutation(double mutationRate = 0.1) : m_mutationRate(mutationRate) {
        if (mutationRate < 0.0 || mutationRate > 1.0) {
            throw std::invalid_argument("Mutation rate must be between 0 and 1");
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

    void setMutationRate(double rate) override {
        m_mutationRate = std::clamp(rate, 0.1, 0.4);
    }

    double getMutationRate() const override {
        return m_mutationRate;
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

    void setMutationRate(double rate) override {
        m_mutationRate = std::clamp(rate, 0.1, 0.4);
    }

    double getMutationRate() const override {
        return m_mutationRate;
    }

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
    double m_mutationRate;
    std::shared_ptr<IMutation> m_pointMutation;
    std::shared_ptr<IMutation> m_swapMutation;
    std::shared_ptr<IMutation> m_guidedMutation;
    
    // Stagnation tracking
    int m_stagnationCounter;
    double m_lastBestFitness;
    double m_adaptiveMutationRate;
    
    // Constants for adaptive mutation
    static constexpr double MIN_MUTATION_RATE = 0.1;
    static constexpr double MAX_MUTATION_RATE = 0.4;
    static constexpr double IMPROVEMENT_THRESHOLD = 0.001;
    static constexpr int STAGNATION_THRESHOLD = 3;

public:
    explicit CombinedMutation(double initialMutationRate = 0.1) 
        : m_mutationRate(initialMutationRate)
        , m_pointMutation(std::make_shared<PointMutation>(initialMutationRate))
        , m_swapMutation(std::make_shared<SwapMutation>(initialMutationRate))
        , m_guidedMutation(std::make_shared<GuidedMutation>(initialMutationRate))
        , m_stagnationCounter(0)
        , m_lastBestFitness(0.0)
        , m_adaptiveMutationRate(initialMutationRate)
    {
        if (initialMutationRate < MIN_MUTATION_RATE || initialMutationRate > MAX_MUTATION_RATE) {
            throw std::invalid_argument("Mutation rate must be between " + 
                std::to_string(MIN_MUTATION_RATE) + " and " + 
                std::to_string(MAX_MUTATION_RATE));
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
        
    void setMutationRate(double rate) override {
        m_mutationRate = std::clamp(rate, MIN_MUTATION_RATE, MAX_MUTATION_RATE);
        m_adaptiveMutationRate = m_mutationRate;
        if (m_pointMutation) m_pointMutation->setMutationRate(m_mutationRate);
        if (m_swapMutation) m_swapMutation->setMutationRate(m_mutationRate);
        if (m_guidedMutation) m_guidedMutation->setMutationRate(m_mutationRate);
    }
    
    double getMutationRate() const override {
        return m_mutationRate;
    }

    // Update mutation rate based on fitness improvement
    void updateMutationRate(double currentBestFitness);
}; 