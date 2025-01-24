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
TournamentSelection::TournamentSelection(GAConfig& config, std::shared_ptr<IPopulationCache> cache) 
        : m_config(config), m_tournamentSize(config.getTournamentSize()), m_cache(cache) {
    LOG_INFO("TournamentSelection initialized with tournament size: " + std::to_string(m_tournamentSize));
}

std::vector<std::shared_ptr<std::vector<int>>> 
TournamentSelection::select(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation)
{
    if (population.empty()) {
        return {};
    }

    std::vector<std::shared_ptr<std::vector<int>>> selected;
    selected.reserve(population.size());

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, population.size() - 1);

    // For each tournament
    for (size_t i = 0; i < population.size(); ++i) {
        double bestFitness = -std::numeric_limits<double>::infinity();
        std::shared_ptr<std::vector<int>> winner;

        // Run tournament
        for (int j = 0; j < m_tournamentSize; ++j) {
            auto candidate = population[dis(gen)];
            double candidateFitness = m_cache->getOrCalculateFitness(candidate, instance, fitness, representation);

            if (candidateFitness > bestFitness) {
                bestFitness = candidateFitness;
                winner = candidate;
            }
        }

        if (winner) {
            selected.push_back(winner);
        }
    }

    return selected;
}
