#pragma once

#include "interfaces/i_algorithm.h"
#include "interfaces/i_representation.h"
#include "utils/random.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include <memory>
#include <vector>
#include <limits>
#include "individual.h"
#include "fitness_impl.h"
#include "representation.h"
#include <functional>
#include <mutex>
#include "../interfaces/i_selection.h"
#include "../interfaces/i_crossover.h"
#include "../interfaces/i_mutation.h"
#include "../interfaces/i_replacement.h"
#include "../interfaces/i_fitness.h"
#include "../interfaces/i_stopping.h"
#include "../interfaces/i_population_cache.h"
#include "adaptive_crossover.h"
#include "path_analyzer.h"
#include "preprocessed_edge.h"

class GeneticAlgorithm : public IAlgorithm {
public:
    GeneticAlgorithm(
        std::unique_ptr<IRepresentation> representation,
        GAConfig& config,
        bool debugMode = false);
    ~GeneticAlgorithm() override = default;

    std::string run(const DNAInstance& instance) override;

    void setProcessId(int id) { m_processId = id; }
    double getBestFitness() const { return m_bestFitness; }
    std::string getBestDNA() const { return m_bestDNA; }

private:
    void initializePopulation(int popSize, const DNAInstance& instance);
    void evaluatePopulation(
        const DNAInstance& instance,
        const std::vector<std::shared_ptr<Individual>>& population,
        const std::shared_ptr<IRepresentation>& representation);
    bool updateGlobalBest(
        const std::vector<std::shared_ptr<Individual>>& population,
        const DNAInstance& instance);
    void logGenerationStats(
        const std::vector<std::shared_ptr<Individual>>& population,
        const DNAInstance& instance,
        int generation);
    double calculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance);
    void calculateTheoreticalMaxFitness(const DNAInstance& instance);
    std::string vectorToString(const std::vector<int>& vec);
    
    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(
        const DNAInstance& instance) const;
    int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;
    int calculateLevenshteinDistance(const std::string& s1, const std::string& s2);

    std::vector<std::shared_ptr<Individual>> m_population;
    std::unique_ptr<IRepresentation> m_representation;
    GAConfig& m_config;
    std::unique_ptr<Random> m_random;
    double m_globalBestFit;
    std::shared_ptr<Individual> m_globalBestIndividual;
    double m_bestFitness;
    double m_theoreticalMaxFitness;
    std::vector<int> m_globalBestGenes;
    std::string m_bestDNA;
    bool m_debugMode;
    int m_processId{0};
    int m_bestIndex;
    std::shared_ptr<IFitness> m_fitness{std::make_shared<OptimizedGraphBasedFitness>()};

    bool shouldStop(int generation) const {
        auto stopping = m_config.getStopping();
        return stopping && stopping->shouldStop(generation, m_globalBestFit);
    }
};
