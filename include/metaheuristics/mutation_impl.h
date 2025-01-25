#pragma once

#include "../interfaces/i_mutation.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include <random>
#include <algorithm>
#include <mutex>

namespace {
    // Thread-safe random number generator
    class RandomGenerator {
    public:
        static RandomGenerator& getInstance() {
            static RandomGenerator instance;
            return instance;
        }

        double getRandomReal(double min = 0.0, double max = 1.0) {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::uniform_real_distribution<double> dist(min, max);
            return dist(m_gen);
        }

        int getRandomInt(int min, int max) {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::uniform_int_distribution<int> dist(min, max);
            return dist(m_gen);
        }

    private:
        RandomGenerator() : m_gen(std::random_device{}()) {}
        std::mt19937 m_gen;
        std::mutex m_mutex;
    };
}

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
            auto& rng = RandomGenerator::getInstance();
            auto genes = individual->getGenes();
            bool mutated = false;

            // Store original genes for rollback if needed
            auto originalGenes = genes;

            // Try to mutate each position with probability m_mutationRate
            for (size_t i = 0; i < genes.size(); ++i) {
                if (rng.getRandomReal() < m_mutationRate) {
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
            if (representation->isValid(genes, instance)) {
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