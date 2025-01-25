#pragma once

#include "../interfaces/i_selection.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "../utils/logging.h"
#include <vector>
#include <memory>
#include <random>
#include <mutex>
#include <algorithm>

namespace {
    // Thread-safe random number generator
    class RandomGenerator {
    public:
        static RandomGenerator& getInstance() {
            static RandomGenerator instance;
            return instance;
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

class TournamentSelection : public ISelection {
public:
    TournamentSelection(GAConfig& config, std::shared_ptr<IPopulationCache> cache)
        : m_config(config)
        , m_tournamentSize(config.getTournamentSize())
        , m_cache(cache)
    {
        if (m_tournamentSize <= 0) {
            LOG_WARNING("Invalid tournament size " + std::to_string(m_tournamentSize) + " - setting to 2");
            m_tournamentSize = 2;
        }
    }
    
    std::vector<std::shared_ptr<Individual>> select(
        const std::vector<std::shared_ptr<Individual>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) override {
        if (population.empty()) {
            LOG_WARNING("Empty population provided to selection operator");
            return {};
        }

        if (!fitness || !representation) {
            LOG_ERROR("Null fitness or representation provided to selection operator");
            return {};
        }

        try {
            std::vector<std::shared_ptr<Individual>> selected;
            selected.reserve(m_config.getParentCount());
            auto& rng = RandomGenerator::getInstance();

            // Run tournaments until we have enough parents
            while (selected.size() < m_config.getParentCount()) {
                // Select tournament participants
                std::vector<std::shared_ptr<Individual>> tournament;
                tournament.reserve(m_tournamentSize);

                for (int i = 0; i < m_tournamentSize; ++i) {
                    int idx = rng.getRandomInt(0, population.size() - 1);
                    if (population[idx] && population[idx]->isValid()) {
                        tournament.push_back(population[idx]);
                    }
                }

                if (tournament.empty()) {
                    LOG_WARNING("No valid individuals in tournament");
                    continue;
                }

                // Find the best individual in the tournament
                auto best = tournament[0];
                double bestFitness = fitness->calculateFitness(best, instance, representation);

                for (size_t i = 1; i < tournament.size(); ++i) {
                    double currentFitness = fitness->calculateFitness(tournament[i], instance, representation);
                    if (currentFitness > bestFitness) {
                        best = tournament[i];
                        bestFitness = currentFitness;
                    }
                }

                // Add the winner to selected parents
                selected.push_back(std::make_shared<Individual>(*best));
            }

            return selected;
        } catch (const std::exception& e) {
            LOG_ERROR("Error during selection: " + std::string(e.what()));
            return {};
        }
    }

private:
    GAConfig& m_config;
    int m_tournamentSize;
    std::shared_ptr<IPopulationCache> m_cache;
    mutable std::mutex m_mutex;  // For thread safety
}; 