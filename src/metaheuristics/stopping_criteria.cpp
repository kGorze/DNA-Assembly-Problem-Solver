//
// Created by konrad_guest on 07/01/2025.
// SMART
#include <chrono>
#include "metaheuristics/stopping_criteria.h"
#include "configuration/genetic_algorithm_configuration.h"
#include <iostream>

bool NoImprovementStopping::stop(const std::vector<std::shared_ptr<std::vector<int>>> &population,
                                 int generation,
                                 const DNAInstance &instance,
                                 std::shared_ptr<IFitness> fitness,
                                 std::shared_ptr<IRepresentation> representation)
{
    double currentBestFitness = -std::numeric_limits<double>::infinity();
    for (auto& individual : population) {
        double fit = fitness->evaluate(individual, instance, representation);
        if (fit > currentBestFitness) {
            currentBestFitness = fit;
        }
    }
    
    if (currentBestFitness > m_bestFitness) {
        m_bestFitness = currentBestFitness;
        m_generationsWithoutImprovement = 0;
    } else {
        m_generationsWithoutImprovement++;
    }
    
    return m_generationsWithoutImprovement >= m_maxGenerationsWithoutImprovement;
}

MaxGenerationsStopping::MaxGenerationsStopping(GAConfig& config) 
    : m_maxGenerations(config.getMaxGenerations())  // Initialize with config value
    , m_useConfig(true)
{
    std::cout << "[MaxGenerationsStopping] Created with maxGenerations = " << m_maxGenerations << " from config" << std::endl;
}

MaxGenerationsStopping::MaxGenerationsStopping(int maxGen)
    : m_maxGenerations(maxGen)
    , m_useConfig(false)
{
    std::cout << "[MaxGenerationsStopping] Created with fixed maxGenerations = " << maxGen << std::endl;
}

bool MaxGenerationsStopping::stop(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    int generation,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation)
{
    // Always use the local value, which was either set from config or explicitly
    std::cout << "[MaxGenerationsStopping] Generation " << generation << " of " << m_maxGenerations 
              << (m_useConfig ? " (from config)" : " (local value)") << std::endl;
    return generation >= m_maxGenerations;
}
