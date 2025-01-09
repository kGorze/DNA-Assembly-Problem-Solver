//
// Created by konrad_guest on 07/01/2025.
//
#include <chrono>
#include "metaheuristics/stopping_criteria.h"

bool NoImprovementStopping::stop(const std::vector<void*>& population,
                  int generation,
                  const DNAInstance& instance,
                  std::shared_ptr<IFitness> fitness,
                  std::shared_ptr<IRepresentation> representation){
        
    // Find best fitness in current population
    double currentBestFitness = -std::numeric_limits<double>::infinity();
    for (auto& individual : population) {
        double fit = fitness->evaluate(individual, instance, representation);
        currentBestFitness = std::max(currentBestFitness, fit);
    }
        
    // If we found a better solution, reset counter
    if (currentBestFitness > m_bestFitness) {
        m_bestFitness = currentBestFitness;
        m_generationsWithoutImprovement = 0;
    } else {
        m_generationsWithoutImprovement++;
    }
        
    // Stop if we haven't improved for 30% of total generations
    return m_generationsWithoutImprovement >= m_maxGenerationsWithoutImprovement;
}

// bool MaxGenerationsStopping::stop(const std::vector<void*>         &population,
//               int                               generation,
//               const DNAInstance                &instance,
//               std::shared_ptr<IFitness>         fitness,
//               std::shared_ptr<IRepresentation>  representation)
// {
//     return (generation >= m_maxGens);
// }

// ========== NoImprovementStopping ==========
// NoImprovementStopping::NoImprovementStopping(int maxNoImprove)
//     : m_maxNoImprove(maxNoImprove),
//       m_bestSoFar(-1e9),
//       m_noImproveCount(0)
// {}
//
// bool NoImprovementStopping::stop(const std::vector<std::vector<double>> &population,
//                                  int generation,
//                                  const IFitness &fitness)
// {
//     // Szukamy najlepszego fit w populacji
//     double best = -1e9;
//     for (auto &ind : population) {
//         double f = fitness.evaluate(ind);
//         if(f > best) {
//             best = f;
//         }
//     }
//     if(best > m_bestSoFar) {
//         m_bestSoFar = best;
//         m_noImproveCount = 0;
//     } else {
//         m_noImproveCount++;
//     }
//     return (m_noImproveCount >= m_maxNoImprove);
// }
//
// // ========== TimeLimitStopping ==========
// TimeLimitStopping::TimeLimitStopping(double limitSec)
//     : m_limitSec(limitSec)
// {
//     m_startTime = std::chrono::duration<double>(
//         std::chrono::high_resolution_clock::now().time_since_epoch()
//     ).count();
// }
//
// bool TimeLimitStopping::stop(const std::vector<std::vector<double>> &population,
//                              int generation,
//                              const IFitness &fitness)
// {
//     double now = std::chrono::duration<double>(
//         std::chrono::high_resolution_clock::now().time_since_epoch()
//     ).count();
//
//     return ((now - m_startTime) >= m_limitSec);
// }
