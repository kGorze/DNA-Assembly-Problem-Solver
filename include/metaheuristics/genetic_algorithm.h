#pragma once

#include "interfaces/i_algorithm.h"
#include "interfaces/i_representation.h"
#include "utils/random.h"
#include "dna/dna_instance.h"
#include <memory>
#include <vector>

struct GeneticConfig {
    int populationSize = 100;
    int maxGenerations = 1000;
    double mutationProbability = 0.1;
    double crossoverProbability = 0.8;
    double targetFitness = 0.95;
    int tournamentSize = 3;
};

class GeneticAlgorithm : public IAlgorithm {
public:
    explicit GeneticAlgorithm(std::unique_ptr<IRepresentation> representation, const GeneticConfig& config = GeneticConfig());
    ~GeneticAlgorithm() override = default;

    std::string run(const DNAInstance& instance) override;

    void setProcessId(int id) { m_processId = id; }
    double getBestFitness() const { return m_bestFitness; }
    std::string getBestDNA() const { return m_bestDNA; }

private:
    void initializePopulation(int popSize, const DNAInstance& instance);
    void evaluatePopulation(const DNAInstance& instance, const std::vector<std::shared_ptr<Individual>>& population);
    void updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance);
    void logGenerationStats(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance, int generation);
    double calculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance);
    void calculateTheoreticalMaxFitness(const DNAInstance& instance);
    std::string vectorToString(const std::vector<int>& vec);
    
    struct PreprocessedEdge {
        int to;
        int weight;
    };
    
    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(const DNAInstance& instance) const;
    int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;

    std::unique_ptr<IRepresentation> m_representation;
    std::unique_ptr<Random> m_random;
    GeneticConfig m_config;
    std::vector<std::shared_ptr<Individual>> m_population;
    std::vector<int> m_globalBestGenes;
    double m_globalBestFit = 0.0;
    int m_processId{0};
    double m_bestFitness{0.0};
    std::string m_bestDNA;
};
