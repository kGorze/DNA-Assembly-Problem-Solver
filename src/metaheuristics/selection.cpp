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
TournamentSelection::select(const std::vector<void*> &population,
                            const DNAInstance &instance,
                            std::shared_ptr<IFitness> fitness,
                            std::shared_ptr<IRepresentation> representation)
{
    // Let's produce as many parents as population.size().
    // We'll do 2 picks at a time = 2 tournaments.
    std::vector<void*> parents;
    parents.reserve(population.size());

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, (int)population.size()-1);

    // Each iteration picks 2 parents, so we do population.size()/2 iterations
    for (size_t i = 0; i < population.size()/2; i++) {
        // Parent 1
        double bestVal1 = -1e9;
        void* bestInd1 = nullptr;
        for(int t=0; t < m_tournamentSize; t++){
            int idx = dist(rng);
            double fv = fitness->evaluate(population[idx], instance, representation);
            if(fv > bestVal1) {
                bestVal1 = fv;
                bestInd1 = population[idx];
            }
        }
        parents.push_back(bestInd1);

        // Parent 2
        double bestVal2 = -1e9;
        void* bestInd2 = nullptr;
        for(int t=0; t < m_tournamentSize; t++){
            int idx = dist(rng);
            double fv = fitness->evaluate(population[idx], instance, representation);
            if(fv > bestVal2) {
                bestVal2 = fv;
                bestInd2 = population[idx];
            }
        }
        parents.push_back(bestInd2);
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
