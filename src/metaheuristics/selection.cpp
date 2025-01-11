//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/selection.h"

#include <algorithm>
#include <random>

// ============== TournamentSelection ==============
TournamentSelection::TournamentSelection(int tournamentSize, std::shared_ptr<IPopulationCache> cache) 
        : m_tournamentSize(tournamentSize), m_cache(cache) {}

std::vector<std::shared_ptr<std::vector<int>>> 
TournamentSelection::select(const std::vector<std::shared_ptr<std::vector<int>>> &population,
                            const DNAInstance &instance,
                            std::shared_ptr<IFitness> fitness,
                            std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> parents;
    parents.reserve(population.size());

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, (int)population.size()-1);

    // Each iteration picks 2 parents, so we do population.size()/2 iterations
    for (size_t i = 0; i < population.size()/2; i++) {
        double bestVal1 = -1e9;
        std::shared_ptr<std::vector<int>> bestInd1 = nullptr;
        for(int t = 0; t < m_tournamentSize; t++){
            int idx = dist(rng);
            double fv = m_cache->getOrCalculateFitness(population[idx], instance, fitness, representation);
            if(fv > bestVal1) {
                bestVal1 = fv;
                bestInd1 = population[idx];
            }
        }
        parents.push_back(bestInd1);

        double bestVal2 = -1e9;
        std::shared_ptr<std::vector<int>> bestInd2 = nullptr;
        for(int t = 0; t < m_tournamentSize; t++){
            int idx = dist(rng);
            double fv = m_cache->getOrCalculateFitness(population[idx], instance, fitness, representation);
            if(fv > bestVal2) {
                bestVal2 = fv;
                bestInd2 = population[idx];
            }
        }
        parents.push_back(bestInd2);
    }
    return parents;
}
