#pragma once

#include "../interfaces/i_stopping.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include <limits>

class NoImprovementStopping : public IStopping {
public:
    explicit NoImprovementStopping(int maxGenerationsWithoutImprovement = 20)
        : m_maxGenerationsWithoutImprovement(maxGenerationsWithoutImprovement)
        , m_bestFitness(std::numeric_limits<double>::lowest())
        , m_generationsWithoutImprovement(0)
    {}

    bool stop(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        int currentGeneration,
        double bestFitness
    ) override;

    void reset();

private:
    int m_maxGenerationsWithoutImprovement;
    double m_bestFitness;
    int m_generationsWithoutImprovement;
};

class MaxGenerationsStopping : public IStopping {
public:
    explicit MaxGenerationsStopping(GAConfig& config);
    explicit MaxGenerationsStopping(int maxGen);

    bool stop(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        int currentGeneration,
        double bestFitness
    ) override;

    void reset();

private:
    int m_maxGenerations;
    bool m_useConfig;
}; 