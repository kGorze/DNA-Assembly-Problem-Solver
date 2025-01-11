//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef SELECTION_H
#define SELECTION_H
#include <vector>
#include <memory>
#include "metaheuristics/fitness.h"
#include "metaheuristics/representation.h"
#include "metaheuristics/population_cache.h"

class ISelection {
public:
    virtual ~ISelection() = default;

    virtual std::vector<std::shared_ptr<std::vector<int>>>
    select(const std::vector<std::shared_ptr<std::vector<int>>> &population,
           const DNAInstance &instance,
           std::shared_ptr<IFitness> fitness,
           std::shared_ptr<IRepresentation> representation) = 0;
};

class TournamentSelection : public ISelection {
public:
    TournamentSelection(int tournamentSize, std::shared_ptr<IPopulationCache> cache);
    std::vector<std::shared_ptr<std::vector<int>>>
    select(const std::vector<std::shared_ptr<std::vector<int>>> &population,
           const DNAInstance &instance,
           std::shared_ptr<IFitness> fitness,
           std::shared_ptr<IRepresentation> representation) override;
private:
    int m_tournamentSize;
    std::shared_ptr<IPopulationCache> m_cache;
};

#endif //SELECTION_H
