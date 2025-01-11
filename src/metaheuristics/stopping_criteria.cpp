//
// Created by konrad_guest on 07/01/2025.
// SMART
#include <chrono>
#include "metaheuristics/stopping_criteria.h"

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
