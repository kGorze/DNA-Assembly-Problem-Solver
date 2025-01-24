//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/selection.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "../../include/generator/dna_generator.h"
#include "utils/logging.h"

#include <algorithm>
#include <random>
#include <sstream>

// ============== TournamentSelection ==============
TournamentSelection::TournamentSelection(GAConfig& config, std::shared_ptr<IPopulationCache> cache) 
        : m_config(config), m_tournamentSize(config.tournamentSize), m_cache(cache) {
    LOG_INFO("TournamentSelection initialized with tournament size: " + std::to_string(m_tournamentSize));
}

std::vector<std::shared_ptr<std::vector<int>>> 
TournamentSelection::select(const std::vector<std::shared_ptr<std::vector<int>>>& population,
                          const DNAInstance& instance,
                          std::shared_ptr<IFitness> fitness,
                          std::shared_ptr<IRepresentation> representation)
{
    if (population.empty()) {
        LOG_ERROR("Empty population in selection");
        return {};
    }

    int tournamentSize = m_tournamentSize;
    
    // Ensure minimum tournament size
    if (tournamentSize <= 0) {
        LOG_WARNING("Invalid tournament size (" + std::to_string(tournamentSize) + "), setting to 2");
        tournamentSize = 2;
    }

    DEBUG_LOG("Starting selection process with population size: " + std::to_string(population.size()));
    
    int targetSize = population.size();
    std::vector<std::shared_ptr<std::vector<int>>> selected;
    selected.reserve(targetSize);
    
    static std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, population.size() - 1);
    
    // Select parents
    for (int i = 0; i < targetSize; i++) {
        double bestFitness = -std::numeric_limits<double>::infinity();
        std::shared_ptr<std::vector<int>> bestIndividual = nullptr;
        
        // Run tournament
        for (int j = 0; j < tournamentSize; j++) {
            int idx = dist(rng);
            auto candidate = population[idx];
            if (!candidate) {
                LOG_WARNING("Null candidate encountered in tournament");
                continue;
            }
            
            double candidateFitness = fitness->evaluate(candidate, instance, representation);
            DEBUG_LOG("Tournament candidate " + std::to_string(j) + " fitness: " + std::to_string(candidateFitness));
            
            if (candidateFitness > bestFitness) {
                bestFitness = candidateFitness;
                bestIndividual = std::make_shared<std::vector<int>>(*candidate);
            }
        }
        
        if (bestIndividual) {
            selected.push_back(bestIndividual);
        } else {
            LOG_ERROR("Failed to select individual in tournament " + std::to_string(i));
        }
    }
    
    LOG_INFO("Selection complete. Selected " + std::to_string(selected.size()) + " individuals");
    return selected;
}
