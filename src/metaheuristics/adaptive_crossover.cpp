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
    m_dnaAlignmentCrossover = std::make_shared<DNAAlignmentCrossover>();
    
    crossovers.emplace_back(CrossoverPerformance(m_orderCrossover));
    crossovers.emplace_back(CrossoverPerformance(m_edgeRecombination));
    crossovers.emplace_back(CrossoverPerformance(m_pmxCrossover));
    crossovers.emplace_back(CrossoverPerformance(m_dnaAlignmentCrossover));
    
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
    
    // Calculate improvement magnitude
    double improvementMagnitude = previousBestFitness - bestSeenFitness;
    bool significantImprovement = improved || (improvementMagnitude > EPSILON);
    
    if (significantImprovement) {
        current.successCount++;
        current.recentSuccessCount++;
        
        // Scale success rate by improvement magnitude
        double improvementFactor = 1.0;
        if (improvementMagnitude > EPSILON) {
            improvementFactor = 1.0 + std::min(1.0, improvementMagnitude);
            bestSeenFitness = previousBestFitness;
        }
        
        // Update recent success rate with improvement scaling
        if (current.recentUsageCount > 0) {
            current.recentSuccessRate = (static_cast<double>(current.recentSuccessCount) / 
                                       current.recentUsageCount) * improvementFactor;
        }
    } else {
        // Penalize lack of improvement
        if (current.recentUsageCount > 0) {
            current.recentSuccessRate = (static_cast<double>(current.recentSuccessCount) / 
                                       current.recentUsageCount) * 0.9;  // 10% penalty
        }
    }
    
    // Update overall success rate with inertia
    if (current.usageCount > 0) {
        double historicalRate = static_cast<double>(current.successCount) / current.usageCount;
        current.successRate = (INERTIA * current.successRate) + 
                            ((1.0 - INERTIA) * current.recentSuccessRate);
        
        // Blend with historical rate to prevent over-specialization
        current.successRate = (0.8 * current.successRate) + (0.2 * historicalRate);
    }
    
    // Periodically adjust probabilities
    if (++generationCount % ADAPTATION_INTERVAL == 0) {
        adjustProbabilities();
    }
}

void AdaptiveCrossover::adjustProbabilities() {
    // Calculate total success rate
    double totalRate = 0.0;
    int activeOperators = 0;
    
    for (auto& c : crossovers) {
        if (c.recentUsageCount >= MIN_TRIALS) {
            totalRate += c.successRate;
            activeOperators++;
        }
    }
    
    if (totalRate < EPSILON || activeOperators == 0) {
        // Reset probabilities if no operator is performing well
        double evenRate = 1.0 / crossovers.size();
        for (auto& c : crossovers) {
            c.successRate = evenRate;
        }
        return;
    }
    
    // Adjust probabilities based on success rates
    for (auto& c : crossovers) {
        if (c.recentUsageCount >= MIN_TRIALS) {
            // Scale probability by success rate
            c.successRate = (c.successRate / totalRate) * 0.9;  // Reserve 10% for exploration
        } else {
            // Give unexplored operators a chance
            c.successRate = 0.1 / (crossovers.size() - activeOperators);
        }
        
        // Ensure minimum probability
        c.successRate = std::max(MIN_PROB, c.successRate);
    }
    
    // Normalize probabilities
    totalRate = 0.0;
    for (auto& c : crossovers) {
        totalRate += c.successRate;
    }
    
    for (auto& c : crossovers) {
        c.successRate /= totalRate;
        
        // Reset recent statistics
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

