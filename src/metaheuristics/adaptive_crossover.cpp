//
// Created by konrad_guest on 10/01/2025.
// SMART
#include "../../include/metaheuristics/adaptive_crossover.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/interfaces/i_crossover.h"
#include "../../include/interfaces/i_representation.h"
#include "../../include/utils/logging.h"
#include "../../include/utils/random.h"
#include <random>
#include <iostream>
#include <algorithm>
#include <limits>

namespace {
    [[maybe_unused]]
    static bool isValidPermutation(const std::vector<int>& vec, size_t size) {
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

AdaptiveCrossover::AdaptiveCrossover(const GAConfig& config)
    : previousBestFitness(-std::numeric_limits<double>::infinity())
    , bestSeenFitness(-std::numeric_limits<double>::infinity())
    , currentCrossoverIndex(0)
    , generationCount(0)
    , INERTIA(0.5)
    , ADAPTATION_INTERVAL(10)
    , MIN_TRIALS(5)
    , MIN_PROB(0.1)
    , m_config(config)
    , m_random(std::make_unique<Random>()) {
    
    // Store crossover operators in member variables to ensure they live as long as AdaptiveCrossover
    m_orderCrossover = std::make_shared<OrderCrossover>();
    m_edgeRecombination = std::make_shared<EdgeRecombination>();
    m_pmxCrossover = std::make_shared<PMXCrossover>();
    
    crossovers.emplace_back(CrossoverPerformance(m_orderCrossover));
    crossovers.emplace_back(CrossoverPerformance(m_edgeRecombination));
    crossovers.emplace_back(CrossoverPerformance(m_pmxCrossover));
    
    // Initialize metrics
    metrics.convergenceGeneration = -1;
    metrics.bestFitness = -std::numeric_limits<double>::infinity();
    metrics.avgFitness = 0.0;
    
    // Initialize success rates evenly
    double initialRate = 1.0 / crossovers.size();
    for (auto& crossover : crossovers) {
        crossover.successRate = initialRate;
    }
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

    if (!hasEnoughData) {
        for (auto& c : crossovers) {
            c.successRate = 1.0 / crossovers.size();
        }
    }
    else {
        double avgSuccessRate = 0.0;
        for (const auto& c : crossovers) {
            avgSuccessRate += c.recentSuccessRate;
        }
        avgSuccessRate /= crossovers.size();

        for (auto& c : crossovers) {
            double newRate = (c.recentSuccessRate / (avgSuccessRate + EPSILON));
            c.successRate = INERTIA * c.successRate + (1.0 - INERTIA) * newRate;
        }
    }

    double totalRate = 0.0;
    for (auto& c : crossovers) {
        totalRate += c.successRate;
    }
    if (totalRate < 1e-9) {
        totalRate = 1.0;
    }

    for (auto& c : crossovers) {
        c.successRate /= totalRate;
        c.successRate = std::max(MIN_PROB, c.successRate);
    }

    double sumRates = 0.0;
    for (auto& c : crossovers) {
        sumRates += c.successRate;
    }
    if (sumRates < 1e-9) {
        sumRates = 1.0;
    }
    for (auto& c : crossovers) {
        c.successRate /= sumRates;
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
    
    currentCrossoverIndex = static_cast<int>(crossovers.size() - 1);
    return crossovers.back().crossover;
}

std::vector<std::shared_ptr<Individual>> AdaptiveCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !representation) {
        LOG_WARNING("Invalid input for adaptive crossover");
        return parents;
    }

    // Select crossover operator based on performance
    auto crossoverOp = selectCrossover();
    if (!crossoverOp) {
        LOG_WARNING("No crossover operator selected");
        return parents;
    }

    // Perform crossover - always keep offspring regardless of validation
    auto offspring = crossoverOp->crossover(parents, instance, representation);
    
    // If crossover produced no offspring (technical error), return parents
    if (offspring.empty()) {
        LOG_WARNING("Crossover produced no offspring - technical error");
        return parents;
    }
    
    // Keep offspring regardless of validation status
    // Let fitness function handle penalties for mismatches
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
