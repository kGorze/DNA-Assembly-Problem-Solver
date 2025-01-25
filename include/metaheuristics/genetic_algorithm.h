#pragma once

#include "interfaces/i_algorithm.h"
#include "interfaces/i_representation.h"
#include "metaheuristics/individual.h"
#include "utils/random.h"
#include <memory>
#include <vector>
#include <string>

struct GeneticConfig {
    int populationSize = 100;
    int maxGenerations = 1000;
    double mutationProbability = 0.1;
    double crossoverProbability = 0.8;
    int tournamentSize = 3;
    double targetFitness = 0.95;
};

class GeneticAlgorithm : public IAlgorithm {
public:
    GeneticAlgorithm(std::unique_ptr<IRepresentation> representation, const GeneticConfig& config);

    ~GeneticAlgorithm() override = default;

    std::string run(const DNAInstance& instance) override;

private:
    std::vector<std::shared_ptr<Individual>> evaluatePopulation(const DNAInstance& instance, const std::vector<std::shared_ptr<Individual>>& population);
    std::vector<std::shared_ptr<Individual>> selectParents(const std::vector<std::shared_ptr<Individual>>& population);
    std::vector<std::shared_ptr<Individual>> performCrossover(const std::vector<std::shared_ptr<Individual>>& parents);
    void mutatePopulation(std::vector<std::shared_ptr<Individual>>& population);
    void updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& population);
    void logGenerationStats(int generation, const std::vector<std::shared_ptr<Individual>>& population);

    std::unique_ptr<IRepresentation> m_representation;
    std::unique_ptr<Random> m_random;
    GeneticConfig m_config;
    double m_globalBestFit;
    std::vector<int> m_globalBestGenes;
};
