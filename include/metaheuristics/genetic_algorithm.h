#pragma once

#include "interfaces/i_algorithm.h"
#include "metaheuristics/i_representation.h"
#include "utils/random.h"
#include <memory>
#include <vector>
#include <string>

struct GeneticConfig {
    int populationSize = 100;
    int maxGenerations = 1000;
    double mutationProbability = 0.1;
    double crossoverProbability = 0.8;
    double targetFitness = 1.0;
};

class GeneticAlgorithm : public IAlgorithm {
public:
    explicit GeneticAlgorithm(std::unique_ptr<IRepresentation> representation, GeneticConfig config = GeneticConfig())
        : m_representation(std::move(representation))
        , m_config(std::move(config))
        , m_random(std::make_unique<Random>())
        , m_globalBestFit(-std::numeric_limits<double>::infinity()) {}

    ~GeneticAlgorithm() override = default;

    std::string run(const DNAInstance& instance) override;

private:
    void evaluatePopulation(std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance);
    void updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance);
    void logGenerationStats(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance, int generation);
    double calculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance);

    std::unique_ptr<IRepresentation> m_representation;
    GeneticConfig m_config;
    std::unique_ptr<Random> m_random;
    double m_globalBestFit;
    std::shared_ptr<Individual> m_globalBestSolution;
};
