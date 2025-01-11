//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H


#include <memory>
#include <vector>
#include <mutex>
#include <omp.h>
#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"
#include "metaheuristics/representation.h"
#include "utils/performance_profilling_framework.h"
#include "metaheuristics/population_cache.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "metaheuristics/adaptive_crossover.h"

/**
 * Klasa integrująca wszystkie komponenty algorytmu ewolucyjnego.
 */
class GeneticAlgorithm {
public:
    GeneticAlgorithm(std::shared_ptr<IRepresentation> representation,
                     std::shared_ptr<ISelection> selection,
                     std::shared_ptr<ICrossover> crossover,
                     std::shared_ptr<IMutation> mutation,
                     std::shared_ptr<IReplacement> replacement,
                     std::shared_ptr<IFitness> fitness,
                     std::shared_ptr<IStopping> stopping,
                     std::shared_ptr<IPopulationCache> cache);

    ~GeneticAlgorithm();

    void run(const DNAInstance &instance);
    std::string getBestDNA() const { return m_bestDNA; }

private:
    void initializePopulation(int popSize, const DNAInstance &instance);
    void logGenerationStats(const std::vector<std::shared_ptr<std::vector<int>>> &pop,
                            const DNAInstance &instance,
                            int generation);

    void updateGlobalBest(const std::vector<std::shared_ptr<std::vector<int>>> &pop,
                          const DNAInstance &instance);
    void calculateTheoreticalMaxFitness(const DNAInstance &instance);
    
    std::vector<std::shared_ptr<std::vector<int>>> population;
    std::shared_ptr<IRepresentation> m_representation;
    std::shared_ptr<ISelection> m_selection;
    std::shared_ptr<ICrossover> m_crossover;
    std::shared_ptr<IMutation> m_mutation;
    std::shared_ptr<IReplacement> m_replacement;
    std::shared_ptr<IFitness> m_fitness;
    std::shared_ptr<IStopping> m_stopping;
    std::shared_ptr<IPopulationCache> m_cache;

    std::string m_bestDNA;
    double m_theoreticalMaxFitness = 0.0;
    double m_globalBestFit = -std::numeric_limits<double>::infinity();
    std::shared_ptr<std::vector<int>> m_globalBestInd;
};

#endif //GENETIC_ALGORITHM_H
