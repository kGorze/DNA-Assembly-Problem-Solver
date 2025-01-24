#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H

#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <string>
#include <limits>
#include <iomanip>
#include <sstream>

#include "metaheuristics/representation.h"
#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"
#include "metaheuristics/population_cache.h"
#include "utils/performance_profilling_framework.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "naive/naive_reconstruction.h"
#include "metaheuristics/adaptive_crossover.h"

void runGeneticAlgorithm(const DNAInstance& instance,
                         const std::string& outputFile = "",
                         int processId = 0,
                         const std::string& difficulty = "Unknown");

double runGeneticAlgorithmWrapper(const DNAInstance& instance);


// Definicja typu funkcji-callbacku
// Teraz z dodatkowymi parametrami coverage, edgeScore, theoreticalMax
using ProgressCallback = std::function<void(int generation,
                                            int maxGenerations,
                                            double bestFitness,
                                            double coverage,
                                            double edgeScore,
                                            double theoreticalMax)>;

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

    void setProgressCallback(ProgressCallback callback) { 
        progressCallback = callback; 
    }
    void setProcessId(int pid) { m_processId = pid; }

    // Po uruchomieniu GA można pobrać najlepsze DNA
    std::string getBestDNA() const { return m_bestDNA; }

private:
    void logGenerationStats(const std::vector<std::shared_ptr<std::vector<int>>>& pop,
                            const DNAInstance& instance,
                            int generation);
    void initializePopulation(int popSize, const DNAInstance &instance);
    void updateGlobalBest(const std::vector<std::shared_ptr<std::vector<int>>> &pop,
                          const DNAInstance &instance);
    void calculateTheoreticalMaxFitness(const DNAInstance &instance);
    void evolve(const DNAInstance& instance);

    double getGlobalBestFitness() const {
        return m_globalBestFit;
    }

private:
    static std::mutex outputMutex;

    std::shared_ptr<IRepresentation> m_representation;
    std::shared_ptr<ISelection> m_selection;
    std::shared_ptr<ICrossover> m_crossover;
    std::shared_ptr<IMutation> m_mutation;
    std::shared_ptr<IReplacement> m_replacement;
    std::shared_ptr<IFitness> m_fitness;
    std::shared_ptr<IStopping> m_stopping;
    std::shared_ptr<IPopulationCache> m_cache;

    std::vector<std::shared_ptr<std::vector<int>>> population;
    
    ProgressCallback progressCallback = nullptr;

    // Najlepszy osobnik
    std::shared_ptr<std::vector<int>> m_globalBestInd;
    double m_globalBestFit;
    std::string m_bestDNA;

    int m_processId = 0;
    double m_theoreticalMaxFitness = 0.0;

};

#endif // GENETIC_ALGORITHM_H
