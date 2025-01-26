//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "../include/metaheuristics/selection_impl.h"
#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../../include/generator/dna_generator.h"
#include "../include/utils/logging.h"

#include <algorithm>
#include <random>
#include <sstream>

// ============== TournamentSelection ==============
TournamentSelection::TournamentSelection(const GAConfig& config, std::shared_ptr<IPopulationCache> cache)
    : m_config(config), m_cache(cache) {
    LOG_INFO("TournamentSelection initialized with tournament size: " + std::to_string(config.getTournamentSize()));
}

std::vector<std::shared_ptr<Individual>> TournamentSelection::select(
    const std::vector<std::shared_ptr<Individual>>& population,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation)
{
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

        std::random_device rd;
        std::mt19937 gen(rd());
        
        while (selected.size() < static_cast<size_t>(m_config.getParentCount())) {
            std::vector<std::shared_ptr<Individual>> tournament;
            tournament.reserve(m_config.getTournamentSize());

            for (int i = 0; i < m_config.getTournamentSize(); ++i) {
                std::uniform_int_distribution<size_t> dist(0, population.size() - 1);
                size_t idx = dist(gen);
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
