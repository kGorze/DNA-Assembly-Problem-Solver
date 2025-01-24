#pragma once

#include "../interfaces/i_stopping.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include <vector>
#include <memory>
#include <limits>

class NoImprovementStopping : public IStopping {
public:
    explicit NoImprovementStopping(int maxGenerationsWithoutImprovement);

    bool stop(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        int currentGeneration,
        double bestFitness
    ) const override;

    void reset();

private:
    int m_maxGenerationsWithoutImprovement;
    mutable int m_generationsWithoutImprovement;
    mutable double m_bestFitness;
};

class MaxGenerationsStopping : public IStopping {
public:
    explicit MaxGenerationsStopping(int maxGenerations);
    explicit MaxGenerationsStopping(GAConfig& config);

    bool stop(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        int currentGeneration,
        double bestFitness
    ) const override;

    void reset() {}

private:
    int m_maxGenerations;
    bool m_useConfig;
}; 