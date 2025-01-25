#pragma once

#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <string>
#include <limits>
#include <random>
#include <sstream>

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
#include "../dna/dna_instance.h"
#include "../metaheuristics/individual.h"
#include "../config/genetic_config.h"

// Progress callback type
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

    // Prevent copying
    GeneticAlgorithm(const GeneticAlgorithm&) = delete;
    GeneticAlgorithm& operator=(const GeneticAlgorithm&) = delete;

    // Allow moving
    GeneticAlgorithm(GeneticAlgorithm&&) noexcept = default;
    GeneticAlgorithm& operator=(GeneticAlgorithm&&) noexcept = default;

    ~GeneticAlgorithm();

    void run(const DNAInstance& instance);

    void setProgressCallback(ProgressCallback callback) { 
        std::lock_guard<std::mutex> lock(m_mutex);
        progressCallback = std::move(callback); 
    }

    void setProcessId(int pid) { 
        std::lock_guard<std::mutex> lock(m_mutex);
        m_processId = std::to_string(pid); 
    }

    std::string getBestDNA() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_bestDNA; 
    }

    double getBestFitness() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_globalBestFit; 
    }

    std::shared_ptr<Individual> getBestIndividual() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_globalBestInd ? std::make_shared<Individual>(*m_globalBestInd) : nullptr;
    }

    void calculateTheoreticalMaxFitness(const DNAInstance& instance);

    static std::string vectorToString(const std::vector<int>& vec) {
        std::stringstream ss;
        for (int val : vec) {
            ss << val << " ";
        }
        return ss.str();
    }

private:
    void logGenerationStats(const std::vector<std::shared_ptr<Individual>>& pop,
                          const DNAInstance& instance,
                          int generation);
    void initializePopulation(int popSize, const DNAInstance& instance);
    void updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& pop,
                         const DNAInstance& instance);
    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(
        const DNAInstance& instance) const;
    int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;

    // Thread safety
    mutable std::mutex m_mutex;
    static std::mutex s_outputMutex;

    // Components (const to prevent modification after construction)
    const std::shared_ptr<IRepresentation> m_representation;
    const std::shared_ptr<ISelection> m_selection;
    const std::shared_ptr<ICrossover> m_crossover;
    const std::shared_ptr<IMutation> m_mutation;
    const std::shared_ptr<IReplacement> m_replacement;
    const std::shared_ptr<IFitness> m_fitness;
    const std::shared_ptr<IPopulationCache> m_cache;
    const std::shared_ptr<IStopping> m_stopping;
    const GeneticConfig& m_config;

    // State
    std::vector<std::shared_ptr<Individual>> m_population;
    mutable double m_globalBestFit{-std::numeric_limits<double>::infinity()};
    mutable std::shared_ptr<Individual> m_globalBestInd{nullptr};
    std::string m_bestDNA;
    std::string m_processId;
    double m_theoreticalMaxFitness{0.0};

    // Progress callback
    ProgressCallback progressCallback;
};
