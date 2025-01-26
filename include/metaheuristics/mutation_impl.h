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
    explicit PointMutation(double mutationRate);

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

    double getMutationRate() const { return m_mutationRate; }
};

class SwapMutation : public IMutation {
private:
    double m_mutationRate;

public:
    explicit SwapMutation(double mutationRate) : m_mutationRate(mutationRate) {}

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        if (!individual || !representation) {
            LOG_WARNING("Null individual or representation provided to mutation operator");
            return;
        }

        auto& rng = Random::instance();
        if (rng.generateProbability() < m_mutationRate) {
            auto genes = individual->getGenes();
            int idx1 = rng.getRandomInt(0, static_cast<int>(genes.size() - 1));
            int idx2 = rng.getRandomInt(0, static_cast<int>(genes.size() - 1));
            std::swap(genes[idx1], genes[idx2]);
            
            // Create a new individual with the mutated genes
            auto mutatedIndividual = std::make_shared<Individual>(genes);
            if (representation->isValid(mutatedIndividual, instance)) {
                individual = mutatedIndividual;
            }
        }
    }
}; 