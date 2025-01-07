//
// Created by konrad_guest on 07/01/2025.
//

#include "metaheuristics/selection.h"

#include <algorithm>
#include <random>

// ============== TournamentSelection ==============
TournamentSelection::TournamentSelection(int tournamentSize)
    : m_tournamentSize(tournamentSize)
{}

std::vector<void*> 
TournamentSelection::select(const std::vector<void*>        &population,
                            const DNAInstance               &instance,
                            std::shared_ptr<IFitness>        fitness,
                            std::shared_ptr<IRepresentation> representation)
{
    // We will produce "parents" the same size as population.
    std::vector<void*> parents;
    parents.reserve(population.size());

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, (int)population.size() - 1);

    for(size_t i=0; i<population.size(); i++){
        double bestVal = -1e9;
        void* bestInd  = nullptr;
        // pick m_tournamentSize random individuals
        for(int t=0; t<m_tournamentSize; t++){
            int idx = dist(rng);
            double fitVal = fitness->evaluate(population[idx], instance, representation);
            if(fitVal > bestVal) {
                bestVal = fitVal;
                bestInd = population[idx];
            }
        }
        // bestInd becomes a parent
        parents.push_back(bestInd);
    }
    return parents;
}

// // ============== RouletteSelection ==============
// std::vector<std::vector<double>> 
// RouletteSelection::select(const std::vector<std::vector<double>> &population,
//                           const IFitness &fitness)
// {
//     // ...
//     return {};
// }
//
// // ============== RankingSelection ==============
// std::vector<std::vector<double>> 
// RankingSelection::select(const std::vector<std::vector<double>> &population,
//                          const IFitness &fitness)
// {
//     // ...
//     return {};
// }
//
// // ============== ElitistSelection ==============
// ElitistSelection::ElitistSelection(int eliteCount)
//     : m_eliteCount(eliteCount)
// {}
//
// std::vector<std::vector<double>> 
// ElitistSelection::select(const std::vector<std::vector<double>> &population,
//                          const IFitness &fitness)
// {
//     // ...
//     return {};
// }
