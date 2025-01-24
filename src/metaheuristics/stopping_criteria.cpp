//
// Created by konrad_guest on 07/01/2025.
// SMART
#include <chrono>
#include "metaheuristics/stopping_criteria.h"
#include "configuration/genetic_algorithm_configuration.h"
#include <iostream>

bool NoImprovementStopping::stop(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    int currentGeneration,
    double bestFitness
) {
    if (bestFitness > m_bestFitness) {
        m_bestFitness = bestFitness;
        m_generationsWithoutImprovement = 0;
        return false;
    }
    
    m_generationsWithoutImprovement++;
    return m_generationsWithoutImprovement >= m_maxGenerationsWithoutImprovement;
}

void NoImprovementStopping::reset() {
    m_bestFitness = std::numeric_limits<double>::lowest();
    m_generationsWithoutImprovement = 0;
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
    const DNAInstance& instance,
    int currentGeneration,
    double bestFitness
) {
    return currentGeneration >= m_maxGenerations;
}

void MaxGenerationsStopping::reset() {
    // Nothing to reset
}
