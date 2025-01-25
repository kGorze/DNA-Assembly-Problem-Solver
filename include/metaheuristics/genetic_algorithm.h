#pragma once

#include "interfaces/i_algorithm.h"
#include "interfaces/i_representation.h"
#include "utils/random.h"
#include <memory>
#include <vector>
#include <mutex>

struct GeneticConfig {
    size_t populationSize = 100;
    size_t maxGenerations = 1000;
    double mutationProbability = 0.1;
    double crossoverProbability = 0.8;
    double targetFitness = 1.0;
};

class GeneticAlgorithm : public IAlgorithm {
public:
    GeneticAlgorithm(const GeneticConfig& config, std::unique_ptr<IRepresentation> representation);
    ~GeneticAlgorithm() override = default;

    std::string run(const DNAInstance& instance) override;

private:
    void evaluatePopulation();
    std::vector<std::shared_ptr<Individual>> selectParents();
    std::vector<std::shared_ptr<Individual>> crossover(const std::vector<std::shared_ptr<Individual>>& parents);
    std::pair<std::shared_ptr<Individual>, std::shared_ptr<Individual>> performCrossover(
        const std::shared_ptr<Individual>& parent1,
        const std::shared_ptr<Individual>& parent2);
    void mutate(std::shared_ptr<Individual>& individual);
    void updateBestSolution();
    double calculateFitness(const std::shared_ptr<Individual>& individual);

    GeneticConfig m_config;
    std::unique_ptr<IRepresentation> m_representation;
    std::unique_ptr<Random> m_random;
    std::vector<std::shared_ptr<Individual>> m_population;
    std::shared_ptr<Individual> m_bestSolution;
    double m_bestFitness = -std::numeric_limits<double>::infinity();

    static std::mutex s_outputMutex;
};
