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
    const std::vector<double>& fitness,
    size_t numParents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (population.empty() || fitness.empty() || numParents == 0) {
        return {};
    }

    std::vector<std::shared_ptr<std::vector<int>>> selected;
    selected.reserve(numParents);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, population.size() - 1);

    // For each parent we need to select
    for (size_t i = 0; i < numParents; ++i) {
        // Select tournament participants
        std::vector<size_t> tournamentIndices;
        tournamentIndices.reserve(m_tournamentSize);
        
        for (int j = 0; j < m_tournamentSize; ++j) {
            tournamentIndices.push_back(dist(gen));
        }

        // Find the best participant
        size_t bestIdx = tournamentIndices[0];
        double bestFitness = fitness[bestIdx];

        for (size_t j = 1; j < tournamentIndices.size(); ++j) {
            size_t idx = tournamentIndices[j];
            if (fitness[idx] > bestFitness) {
                bestIdx = idx;
                bestFitness = fitness[idx];
            }
        }

        // Add the winner to selected parents
        selected.push_back(population[bestIdx]);
    }

    return selected;
}
