//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "../include/metaheuristics/stopping_criteria_impl.h"
#include "../include/utils/logging.h"
#include <chrono>
#include <mutex>
#include <sstream>
#include "configuration/genetic_algorithm_configuration.h"
#include <iostream>
#include <limits>

NoImprovementStopping::NoImprovementStopping(int maxGenerationsWithoutImprovement)
    : m_maxGenerationsWithoutImprovement(maxGenerationsWithoutImprovement)
    , m_bestFitness(-std::numeric_limits<double>::infinity())
    , m_generationsWithoutImprovement(0) {}

bool NoImprovementStopping::stop(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    int generation,
    double bestFitness) const {
    if (bestFitness > m_bestFitness) {
        m_bestFitness = bestFitness;
        m_generationsWithoutImprovement = 0;
        return false;
    }
    
    m_generationsWithoutImprovement++;
    return m_generationsWithoutImprovement >= m_maxGenerationsWithoutImprovement;
}

void NoImprovementStopping::reset() {
    m_bestFitness = -std::numeric_limits<double>::infinity();
    m_generationsWithoutImprovement = 0;
}

MaxGenerationsStopping::MaxGenerationsStopping(GAConfig& config) 
    : m_maxGenerations(config.getMaxGenerations())
    , m_useConfig(true)
{
    if (m_maxGenerations <= 0) {
        LOG_WARNING("Invalid maxGenerations value: " + std::to_string(m_maxGenerations) + ". Using default value of 100");
        m_maxGenerations = 100;
    }
    LOG_INFO("MaxGenerationsStopping initialized with maxGenerations = " + std::to_string(m_maxGenerations));
}

MaxGenerationsStopping::MaxGenerationsStopping(int maxGen)
    : m_maxGenerations(maxGen)
    , m_useConfig(false)
{
    if (m_maxGenerations <= 0) {
        LOG_WARNING("Invalid maxGenerations value: " + std::to_string(maxGen) + ". Using default value of 100");
        m_maxGenerations = 100;
    }
    LOG_INFO("MaxGenerationsStopping initialized with fixed maxGenerations = " + std::to_string(m_maxGenerations));
}

bool MaxGenerationsStopping::stop(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    int currentGeneration,
    double bestFitness
) const {
    return currentGeneration >= m_maxGenerations;
}
