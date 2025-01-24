#pragma once

// Standard library includes
#include <memory>
#include <vector>
#include <functional>
#include <mutex>

// Project interface includes
#include "interfaces/i_representation.h"
#include "interfaces/i_selection.h"
#include "interfaces/i_crossover.h"
#include "interfaces/i_mutation.h"
#include "interfaces/i_replacement.h"
#include "interfaces/i_fitness.h"
#include "interfaces/i_stopping.h"
#include "interfaces/i_population_cache.h"

// Project implementation includes
#include "configuration/genetic_algorithm_configuration.h"
#include "dna/dna_instance.h"
#include "metaheuristics/adaptive_crossover.h"
#include "metaheuristics/path_analyzer.h"

class GeneticAlgorithm {
public:
    GeneticAlgorithm(
        std::shared_ptr<IRepresentation> representation,
        std::shared_ptr<ISelection> selection,
        std::shared_ptr<ICrossover> crossover,
        std::shared_ptr<IMutation> mutation,
        std::shared_ptr<IReplacement> replacement,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IStopping> stopping,
        std::shared_ptr<IPopulationCache> cache,
        GAConfig& config
    );

    void run(const DNAInstance& instance);
    std::shared_ptr<std::vector<int>> getBestIndividual() const;
    std::vector<int> getBestDNA() const;
    double getBestFitness() const;
    void setProgressCallback(const std::function<void(int, int, double, double, double, double)>& callback);

private:
    void updateGlobalBest(const std::vector<std::shared_ptr<std::vector<int>>>& pop, const DNAInstance& instance);

    std::shared_ptr<IRepresentation> m_representation;
    std::shared_ptr<ISelection> m_selection;
    std::shared_ptr<ICrossover> m_crossover;
    std::shared_ptr<IMutation> m_mutation;
    std::shared_ptr<IReplacement> m_replacement;
    std::shared_ptr<IFitness> m_fitness;
    std::shared_ptr<IStopping> m_stopping;
    std::shared_ptr<IPopulationCache> m_cache;
    GAConfig& m_config;

    std::vector<std::shared_ptr<std::vector<int>>> population;
    std::shared_ptr<std::vector<int>> m_globalBestInd;
    std::vector<int> m_bestDNA;
    double m_globalBestFit;
    std::function<void(int, int, double, double, double, double)> progressCallback;

    static std::mutex outputMutex;
}; 