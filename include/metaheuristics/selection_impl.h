#pragma once

#include "../interfaces/i_selection.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "../utils/logging.h"
#include "utils/random.h"
#include <vector>
#include <memory>
#include <random>
#include <mutex>
#include <algorithm>

class TournamentSelection : public ISelection {
private:
    const GAConfig& m_config;

public:
    explicit TournamentSelection(const GAConfig& config) : m_config(config) {}

    std::vector<std::shared_ptr<Individual>> select(
            const std::vector<std::shared_ptr<Individual>>& population,
            const DNAInstance& instance,
            std::shared_ptr<IFitness> fitness,
            std::shared_ptr<IRepresentation> representation) override {
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
            selected.reserve(static_cast<size_t>(m_config.getParentCount()));

            auto& rng = Random::instance();
            while (selected.size() < static_cast<size_t>(m_config.getParentCount())) {
                std::vector<std::shared_ptr<Individual>> tournament;
                tournament.reserve(m_config.getTournamentSize());

                for (int i = 0; i < m_config.getTournamentSize(); ++i) {
                    int idx = rng.getRandomInt(0, static_cast<int>(population.size() - 1));
                    tournament.push_back(population[idx]);
                }

                if (tournament.empty()) {
                    LOG_WARNING("No valid individuals in tournament");
                    continue;
                }

                auto best = tournament[0];
                double bestFitness = fitness->calculateFitness(best, instance, representation);

                for (size_t i = 1; i < tournament.size(); ++i) {
                    double currentFitness = fitness->calculateFitness(tournament[i], instance, representation);
                    if (currentFitness > bestFitness) {
                        best = tournament[i];
                        bestFitness = currentFitness;
                    }
                }

                selected.push_back(best);
            }

            return selected;
        } catch (const std::exception& e) {
            LOG_ERROR("Error during selection: " + std::string(e.what()));
            return {};
        }
    }
}; 