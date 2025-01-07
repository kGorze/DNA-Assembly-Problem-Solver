//
// Created by konrad_guest on 07/01/2025.
//

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#pragma once

#include <memory>
#include <vector>
#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"
#include "metaheuristics/representation.h"

/**
 * Klasa integrująca wszystkie komponenty algorytmu ewolucyjnego.
 */
class GeneticAlgorithm {
public:
    GeneticAlgorithm(std::shared_ptr<IRepresentation> representation,
                     std::shared_ptr<ISelection>    selection,
                     std::shared_ptr<ICrossover>   crossover,
                     std::shared_ptr<IMutation>    mutation,
                     std::shared_ptr<IReplacement> replacement,
                     std::shared_ptr<IFitness>     fitness,
                     std::shared_ptr<IStopping>    stopping);

    void run(const DNAInstance &instance); 

private:
    std::shared_ptr<IRepresentation> m_representation;
    std::shared_ptr<ISelection>      m_selection;
    std::shared_ptr<ICrossover>      m_crossover;
    std::shared_ptr<IMutation>       m_mutation;
    std::shared_ptr<IReplacement>    m_replacement;
    std::shared_ptr<IFitness>        m_fitness;
    std::shared_ptr<IStopping>       m_stopping;

    std::vector<void*> population;

    void initializePopulation(int popSize, const DNAInstance &instance);
};



#endif //GENETIC_ALGORITHM_H
