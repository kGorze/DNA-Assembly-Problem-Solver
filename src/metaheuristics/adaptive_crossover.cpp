//
// Created by konrad_guest on 10/01/2025.
// SMART
#include "../../include/metaheuristics/adaptive_crossover.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/interfaces/i_crossover.h"
#include "../../include/interfaces/i_representation.h"
#include "../../include/utils/logging.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <limits>

namespace {
    bool isValidPermutation(const std::vector<int>& vec, size_t size) {
        if (vec.size() != size) return false;
        std::vector<bool> used(size, false);
        for (int val : vec) {
            if (val < 0 || val >= static_cast<int>(size) || used[val]) return false;
            used[val] = true;
        }
        return true;
    }
}

AdaptiveCrossover::AdaptiveCrossover()
    : previousBestFitness(-std::numeric_limits<double>::infinity())
    , bestSeenFitness(-std::numeric_limits<double>::infinity())
    , currentCrossoverIndex(0)
    , generationCount(0)
    , INERTIA(0.7)
    , ADAPTATION_INTERVAL(20)
    , MIN_TRIALS(5)
    , MIN_PROB(0.1)
{
    crossovers.emplace_back(CrossoverPerformance(std::make_shared<OrderCrossover>()));
    crossovers.emplace_back(CrossoverPerformance(std::make_shared<EdgeRecombination>()));
    crossovers.emplace_back(CrossoverPerformance(std::make_shared<PMXCrossover>()));
    
    metrics.convergenceGeneration = -1;
    metrics.bestFitness = -std::numeric_limits<double>::infinity();
    metrics.avgFitness = 0.0;
}

void AdaptiveCrossover::setParameters(double inertia, int adaptInterval, int minTrials, double minProb) {
    INERTIA = inertia;
    ADAPTATION_INTERVAL = adaptInterval;
    MIN_TRIALS = minTrials;
    MIN_PROB = minProb;
}

RunMetrics AdaptiveCrossover::getMetrics() const {
    RunMetrics result;
    result.avgFitness = metrics.avgFitness;
    result.bestFitness = metrics.bestFitness;
    result.convergenceTime = metrics.convergenceGeneration;
    
    result.operatorUsageRates.resize(crossovers.size());
    result.operatorSuccessRates.resize(crossovers.size());
    
    for (size_t i = 0; i < crossovers.size(); i++) {
        result.operatorUsageRates[i] = (generationCount == 0)
            ? 0.0
            : static_cast<double>(crossovers[i].usageCount) / generationCount;

        if (crossovers[i].usageCount > 0) {
            result.operatorSuccessRates[i] = static_cast<double>(crossovers[i].successCount)
                                           / crossovers[i].usageCount;
        } else {
            result.operatorSuccessRates[i] = 0.0;
        }
    }
    
    return result;
}

void AdaptiveCrossover::updatePerformance(bool improved) {
    auto& current = crossovers[currentCrossoverIndex];
    current.usageCount++;
    current.recentUsageCount++;
    
    bool significantImprovement = improved || 
        ((previousBestFitness - bestSeenFitness) > EPSILON);
    
    if (significantImprovement) {
        current.successCount++;
        current.recentSuccessCount++;
        
        if (previousBestFitness - bestSeenFitness > EPSILON) {
            bestSeenFitness = previousBestFitness;
        }
    }
    
    if (current.recentUsageCount > 0) {
        current.recentSuccessRate = static_cast<double>(current.recentSuccessCount) / 
                                    current.recentUsageCount;
    }
}

void AdaptiveCrossover::adjustProbabilities() {
    bool hasEnoughData = true;
    for (const auto& c : crossovers) {
        if (c.recentUsageCount < MIN_TRIALS) {
            hasEnoughData = false;
            break;
        }
    }

    // Jeśli nie mamy wystarczająco danych, ustaw równo
    if (!hasEnoughData) {
        for (auto& c : crossovers) {
            c.successRate = 1.0 / crossovers.size();
        }
    }
    else {
        // Oblicz średni recentSuccessRate
        double avgSuccessRate = 0.0;
        for (const auto& c : crossovers) {
            avgSuccessRate += c.recentSuccessRate;
        }
        avgSuccessRate /= crossovers.size();

        // Aktualizuj successRate z uwzględnieniem INERTIA
        for (auto& c : crossovers) {
            double newRate = (c.recentSuccessRate / (avgSuccessRate + EPSILON));
            c.successRate = INERTIA * c.successRate + (1.0 - INERTIA) * newRate;
        }
    }

    // Normalizacja z uwzględnieniem MIN_PROB
    double totalRate = 0.0;
    for (auto& c : crossovers) {
        totalRate += c.successRate;
    }
    if (totalRate < 1e-9) {
        totalRate = 1.0;
    }

    // Wymuszamy minimalne prawdopodobieństwo i ponownie normalizujemy
    for (auto& c : crossovers) {
        // Najpierw ratio
        c.successRate /= totalRate;
        // Potem MIN_PROB
        c.successRate = std::max(MIN_PROB, c.successRate);
    }

    // Druga faza normalizacji
    double sumRates = 0.0;
    for (auto& c : crossovers) {
        sumRates += c.successRate;
    }
    if (sumRates < 1e-9) {
        sumRates = 1.0;
    }
    for (auto& c : crossovers) {
        c.successRate /= sumRates;
        // Reset statystyk
        c.recentUsageCount = 0;
        c.recentSuccessCount = 0;
        c.recentSuccessRate = 0.0;
    }
}

std::shared_ptr<ICrossover> AdaptiveCrossover::selectCrossover() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    double randVal = dis(gen);
    double sum = 0.0;
    
    for (size_t i = 0; i < crossovers.size(); i++) {
        sum += crossovers[i].successRate;
        if (randVal <= sum) {
            currentCrossoverIndex = static_cast<int>(i);
            return crossovers[i].crossover;
        }
    }
    
    // Zapobiegawczo, jeśli nie wylosowano wcześniej
    currentCrossoverIndex = static_cast<int>(crossovers.size() - 1);
    return crossovers.back().crossover;
}

std::vector<std::shared_ptr<std::vector<int>>> AdaptiveCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) 
{
    if (parents.empty()) {
        std::cerr << "[ERROR] Empty parents in crossover\n";
        return {};
    }
    
    auto selectedCrossover = selectCrossover();
    if (!selectedCrossover) {
        return {};
    }
    auto offspring = selectedCrossover->crossover(parents, instance, representation);
    int startIndex = instance.getStartIndex();
    
    // Validate and fix offspring
    for (auto it = offspring.begin(); it != offspring.end();) {
        if (!(*it) || (*it)->empty() || !isValidPermutation(**it, instance.getSpectrum().size())) {
            it = offspring.erase(it);
        } else {
            // Ensure startIndex is at the beginning
            auto startIt = std::find((*it)->begin(), (*it)->end(), startIndex);
            if (startIt != (*it)->begin()) {
                std::iter_swap((*it)->begin(), startIt);
            }
            ++it;
        }
    }
    
    generationCount++;
    
    if (generationCount % ADAPTATION_INTERVAL == 0) {
        adjustProbabilities();
    }
    return offspring;
}

void AdaptiveCrossover::updateFeedback(double currentBestFitness) {
    bool improved = (currentBestFitness > previousBestFitness);
    updatePerformance(improved);
    previousBestFitness = currentBestFitness;
    
    metrics.fitnessHistory.push_back(currentBestFitness);
    if (improved && metrics.convergenceGeneration == -1) {
        metrics.convergenceGeneration = generationCount;
    }
    metrics.bestFitness = std::max(metrics.bestFitness, currentBestFitness);
    if (generationCount > 0) {
        metrics.avgFitness = 
          (metrics.avgFitness * (generationCount - 1) + currentBestFitness)
          / generationCount;
    }
}
