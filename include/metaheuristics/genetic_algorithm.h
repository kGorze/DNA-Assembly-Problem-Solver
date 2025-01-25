#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <string>
#include <limits>
#include <iomanip>
#include <sstream>
#include <random>
#include <algorithm>
#include <chrono>

// Interface includes
#include "../interfaces/i_representation.h"
#include "../interfaces/i_selection.h"
#include "../interfaces/i_crossover.h"
#include "../interfaces/i_mutation.h"
#include "../interfaces/i_replacement.h"
#include "../interfaces/i_fitness.h"
#include "../interfaces/i_stopping.h"
#include "../interfaces/i_population_cache.h"

// Other includes
#include "../utils/performance_profilling_framework.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "../naive/naive_reconstruction.h"
#include "../dna/dna_instance.h"
#include "../metaheuristics/adaptive_crossover.h"
#include "../metaheuristics/path_analyzer.h"
#include "../metaheuristics/preprocessed_edge.h"
#include "../metaheuristics/individual.h"
#include "../config/genetic_config.h"

// Forward declarations of interfaces
class IRepresentation;
class ISelection;
class ICrossover;
class IMutation;
class IReplacement;
class IFitness;
class IStopping;
class IPopulationCache;

double runGeneticAlgorithmWrapper(const DNAInstance& instance);

// Definicja typu funkcji-callbacku
// Teraz z dodatkowymi parametrami coverage, edgeScore, theoreticalMax
using ProgressCallback = std::function<void(int, int, double, double, double, double)>;

class GeneticAlgorithm {
public:
    GeneticAlgorithm(
        std::shared_ptr<IRepresentation> representation,
        std::shared_ptr<ISelection> selection,
        std::shared_ptr<ICrossover> crossover,
        std::shared_ptr<IMutation> mutation,
        std::shared_ptr<IReplacement> replacement,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IPopulationCache> cache,
        std::shared_ptr<IStopping> stopping,
        const GeneticConfig& config
    );

    ~GeneticAlgorithm();

    void run(const DNAInstance& instance);

    void setProgressCallback(ProgressCallback callback) { 
        progressCallback = callback; 
    }
    void setProcessId(int pid) { m_processId = pid; }

    // Po uruchomieniu GA można pobrać najlepsze DNA
    const std::string& getBestDNA() const { return m_bestDNA; }

    double getBestFitness() const { return m_globalBestFit; }
    std::shared_ptr<Individual> getBestIndividual() const { return m_globalBestInd; }

private:
    void logGenerationStats(const std::vector<std::shared_ptr<Individual>>& pop,
                          const DNAInstance& instance,
                          int generation);
    void initializePopulation(int popSize, const DNAInstance& instance);
    void updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& pop,
                         const DNAInstance& instance);
    void calculateTheoreticalMaxFitness(const DNAInstance& instance);
    void evolve(const DNAInstance& instance);
    std::string vectorToString(const std::vector<int>& vec);
    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(const DNAInstance& instance) const;
    int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;

    static std::mutex outputMutex;

    std::shared_ptr<IRepresentation> m_representation;
    std::shared_ptr<ISelection> m_selection;
    std::shared_ptr<ICrossover> m_crossover;
    std::shared_ptr<IMutation> m_mutation;
    std::shared_ptr<IReplacement> m_replacement;
    std::shared_ptr<IFitness> m_fitness;
    std::shared_ptr<IPopulationCache> m_cache;
    std::shared_ptr<IStopping> m_stopping;
    GeneticConfig m_config;

    std::vector<std::shared_ptr<Individual>> population;
    
    ProgressCallback progressCallback;

    // Najlepszy osobnik
    std::shared_ptr<Individual> m_globalBestInd;
    double m_globalBestFit = -std::numeric_limits<double>::infinity();
    std::string m_bestDNA;

    int m_processId = 0;
    double m_theoreticalMaxFitness = 0.0;

    // Random number generation
    std::mt19937 generator;
    std::uniform_real_distribution<> distribution;
};
