#pragma once

#include "../interfaces/i_mutation.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include "utils/random.h"
#include <random>
#include <algorithm>
#include <mutex>

class PointMutation : public IMutation {
public:
    explicit PointMutation(double mutationRate) 
        : m_mutationRate(mutationRate)
    {
        if (m_mutationRate < 0.0 || m_mutationRate > 1.0) {
            LOG_WARNING("Invalid mutation rate " + std::to_string(m_mutationRate) + " - clamping to valid range");
            m_mutationRate = std::clamp(m_mutationRate, 0.0, 1.0);
        }
    }

    void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override {
        if (!individual || !representation) {
            LOG_WARNING("Null individual or representation provided to mutation operator");
            return;
        }

        try {
            auto& rng = Random::instance();
            auto genes = individual->getGenes();
            bool mutated = false;

            // Store original genes for rollback if needed
            auto originalGenes = genes;

            // Try to mutate each position with probability m_mutationRate
            for (size_t i = 0; i < genes.size(); ++i) {
                if (rng.generateProbability() < m_mutationRate) {
                    int oldValue = genes[i];
                    int newValue;
                    do {
                        newValue = rng.getRandomInt(0, instance.getSpectrum().size() - 1);
                    } while (newValue == oldValue);
                    
                    genes[i] = newValue;
                    mutated = true;

                    LOG_DEBUG("Mutated position " + std::to_string(i) + 
                             " from " + std::to_string(oldValue) + 
                             " to " + std::to_string(newValue));
                }
            }

            // Force at least one mutation if none occurred
            if (!mutated) {
                size_t pos = rng.getRandomInt(0, genes.size() - 1);
                int oldValue = genes[pos];
                int newValue;
                do {
                    newValue = rng.getRandomInt(0, instance.getSpectrum().size() - 1);
                } while (newValue == oldValue);
                
                genes[pos] = newValue;
                LOG_DEBUG("Forced mutation at position " + std::to_string(pos) + 
                         " from " + std::to_string(oldValue) + 
                         " to " + std::to_string(newValue));
            }

            // Validate and update the individual
            auto mutatedIndividual = std::make_shared<Individual>();
            mutatedIndividual->setGenes(genes);
            if (representation->isValid(mutatedIndividual, instance)) {
                individual->setGenes(std::move(genes));
            } else {
                LOG_WARNING("Mutation produced invalid solution - rolling back changes");
                individual->setGenes(std::move(originalGenes));
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Error during mutation: " + std::string(e.what()));
        }
    }

    double getMutationRate() const { return m_mutationRate; }

private:
    double m_mutationRate;
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
            
            // Validate and update the individual
            auto mutatedIndividual = std::make_shared<Individual>();
            mutatedIndividual->setGenes(genes);
            if (representation->isValid(mutatedIndividual, instance)) {
                individual->setGenes(std::move(genes));
            }
        }
    }
}; 