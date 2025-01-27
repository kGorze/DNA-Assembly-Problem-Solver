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
#include <sstream>
#include <iomanip>

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

AdaptiveCrossover::AdaptiveCrossover(const GAConfig& config, const DNAInstance& instance)
    : m_config(config)
    , m_instance(instance)
    , m_random(std::make_unique<Random>())
    , m_lastDiversityMeasure(0.0)
    , previousBestFitness(-std::numeric_limits<double>::infinity())
    , bestSeenFitness(-std::numeric_limits<double>::infinity())
    , currentCrossoverIndex(0)
    , generationCount(0)
    , INERTIA(0.5)
    , ADAPTATION_INTERVAL(10)
    , MIN_TRIALS(5)
    , MIN_PROB(0.1) {
    
    // Initialize crossovers
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
    
    LOG_DEBUG("AdaptiveCrossover initialized with ADAPTATION_INTERVAL=" + 
              std::to_string(ADAPTATION_INTERVAL));
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
    LOG_DEBUG("Entering updatePerformance, generation: " + std::to_string(generationCount));
    
    // Update performance metrics for current crossover operator
    if (currentCrossoverIndex >= 0 && currentCrossoverIndex < static_cast<int>(crossovers.size())) {
        auto& currentOp = crossovers[currentCrossoverIndex];
        currentOp.trials++;
        currentOp.recentUsageCount++;
        if (improved) {
            currentOp.successes++;
            currentOp.recentSuccessCount++;
        }
        
        // Update success rates
        if (currentOp.trials > 0) {
            currentOp.successRate = static_cast<double>(currentOp.successes) / currentOp.trials;
        }
        if (currentOp.recentUsageCount > 0) {
            currentOp.recentSuccessRate = static_cast<double>(currentOp.recentSuccessCount) / currentOp.recentUsageCount;
        }
        
        // Log operator performance
        LOG_DEBUG("Operator " + std::to_string(currentCrossoverIndex) + 
                 " performance - Success rate: " + std::to_string(currentOp.successRate) + 
                 ", Recent success rate: " + std::to_string(currentOp.recentSuccessRate));
    }
    
    // Log diversity metrics if we have a cache
    if (m_config.getCache()) {
        auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
        if (cache) {
            const auto& population = cache->getCurrentPopulation();
            auto representation = m_config.getRepresentation();
            if (representation) {
                logDiversityMetrics(population, representation, metrics);
            }
        }
    }
    
    // Store operator usage and success rates
    std::vector<double> usageRates;
    std::vector<double> successRates;
    for (const auto& op : crossovers) {
        usageRates.push_back(static_cast<double>(op.recentUsageCount) / std::max(1, ADAPTATION_INTERVAL));
        successRates.push_back(op.recentSuccessRate);
    }
    metrics.operatorUsageHistory.push_back(usageRates);
    metrics.operatorSuccessHistory.push_back(successRates);
    
    // Reset recent counters if adaptation interval is reached
    if (generationCount % ADAPTATION_INTERVAL == 0) {
        for (auto& op : crossovers) {
            op.recentUsageCount = 0;
            op.recentSuccessCount = 0;
            op.recentSuccessRate = 0.0;
        }
        
        // Log diversity metrics after resetting counters
        if (m_config.getCache()) {
            auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
            if (cache) {
                const auto& population = cache->getCurrentPopulation();
                auto representation = m_config.getRepresentation();
                if (representation) {
                    logDiversityMetrics(population, representation, metrics);
                }
            }
        }
        
        // Adjust probabilities based on performance
        adjustProbabilities();
    }
}

void AdaptiveCrossover::adjustProbabilities() {
    LOG_DEBUG("Adjusting probabilities based on performance");
    
    // Calculate total success rate
    double totalSuccessRate = 0.0;
    for (const auto& op : crossovers) {
        totalSuccessRate += op.successRate;
    }
    
    // Adjust probabilities based on success rates
    if (totalSuccessRate > EPSILON) {
        for (auto& op : crossovers) {
            double newRate = (1.0 - MIN_PROB * crossovers.size()) * 
                           (op.successRate / totalSuccessRate) + MIN_PROB;
            op.successRate = INERTIA * op.successRate + (1.0 - INERTIA) * newRate;
        }
        
        // Log diversity metrics after probability adjustment
        if (m_config.getCache()) {
            auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
            if (cache) {
                const auto& population = cache->getCurrentPopulation();
                auto representation = m_config.getRepresentation();
                if (representation) {
                    logDiversityMetrics(population, representation, metrics);
                }
            }
        }
        
        // Normalize probabilities
        double sum = 0.0;
        for (const auto& op : crossovers) {
            sum += op.successRate;
        }
        if (sum > EPSILON) {
            for (auto& op : crossovers) {
                op.successRate /= sum;
            }
        }
    } else {
        // If no success, use uniform distribution
        double uniformProb = 1.0 / crossovers.size();
        for (auto& op : crossovers) {
            op.successRate = uniformProb;
        }
    }
    
    // Log adjusted probabilities
    std::stringstream ss;
    ss << "Adjusted probabilities:";
    for (size_t i = 0; i < crossovers.size(); ++i) {
        ss << "\n  Operator " << i << ": " << crossovers[i].successRate;
    }
    LOG_DEBUG(ss.str());
}

void AdaptiveCrossover::logDiversityMetrics(
    const std::vector<std::shared_ptr<Individual>>& population,
    std::shared_ptr<IRepresentation> representation,
    Metrics& metrics) {
    
    LOG_DEBUG("Attempting to log diversity metrics at generation " + std::to_string(generationCount));
    
    if (population.empty()) {
        LOG_WARNING("Empty population, skipping diversity metrics");
        return;
    }
    
    // Calculate average Hamming distance
    double avgDistance = calculateAverageDistance(population);
    
    // Calculate unique solutions ratio
    std::unordered_set<std::string> uniqueGenes;
    std::unordered_set<std::string> uniqueDNA;
    
    for (const auto& individual : population) {
        if (!individual) {
            LOG_DEBUG("Skipping null individual");
            continue;
        }
        
        // Use representation's validation
        if (!representation->isValid(individual, m_instance)) {
            LOG_DEBUG("Skipping individual with invalid genes");
            continue;
        }
        
        // Add to unique genes set
        uniqueGenes.insert(representation->toString(individual, m_instance));
        
        // Convert to DNA if possible
        if (representation) {
            try {
                std::string dna = representation->toDNA(individual, m_instance);
                if (!dna.empty()) {
                    uniqueDNA.insert(dna);
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error converting individual to DNA: " + std::string(e.what()));
                continue;
            }
        }
    }
    
    double uniqueSolutionsRatio = population.empty() ? 0.0 : 
                                 static_cast<double>(uniqueGenes.size()) / population.size();
    double uniqueDNASequencesRatio = population.empty() ? 0.0 : 
                                    static_cast<double>(uniqueDNA.size()) / population.size();
    
    // Log all metrics with detailed information
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4);
    ss << "Generation " << generationCount << " Diversity Metrics:";
    ss << "\n  Average Hamming Distance: " << avgDistance;
    ss << "\n  Unique Solutions Ratio: " << uniqueSolutionsRatio << " (" << uniqueGenes.size() << "/" << population.size() << ")";
    ss << "\n  Unique DNA Sequences Ratio: " << uniqueDNASequencesRatio << " (" << uniqueDNA.size() << "/" << population.size() << ")";
    
    LOG_INFO(ss.str());
    
    // Store metrics for later use
    m_lastDiversityMeasure = avgDistance;
    DiversityMetrics diversityMetrics = {avgDistance, uniqueSolutionsRatio, uniqueDNASequencesRatio};
    metrics.diversityHistory.push_back(diversityMetrics);
    
    LOG_DEBUG("Successfully logged diversity metrics for generation " + std::to_string(generationCount));
}

double AdaptiveCrossover::calculateAverageDistance(
    const std::vector<std::shared_ptr<Individual>>& population) const {
    if (population.size() < 2) return 0.0;
    
    double totalDistance = 0.0;
    int comparisons = 0;
    
    for (size_t i = 0; i < population.size(); ++i) {
        for (size_t j = i + 1; j < population.size(); ++j) {
            if (!population[i] || !population[j]) continue;
            
            const auto& genes1 = population[i]->getGenes();
            const auto& genes2 = population[j]->getGenes();
            
            // Calculate normalized Hamming distance
            double distance = 0.0;
            size_t minSize = std::min(genes1.size(), genes2.size());
            size_t maxSize = std::max(genes1.size(), genes2.size());
            
            for (size_t k = 0; k < minSize; ++k) {
                if (genes1[k] != genes2[k]) distance += 1.0;
            }
            
            // Add penalty for length difference
            distance += (maxSize - minSize);
            
            // Normalize
            distance /= maxSize;
            
            totalDistance += distance;
            comparisons++;
        }
    }
    
    return comparisons > 0 ? totalDistance / comparisons : 0.0;
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
    
    LOG_DEBUG("Entering crossover, generation: " + std::to_string(generationCount));

    if (parents.size() < 2 || !representation) {
        LOG_ERROR("Invalid parents or representation in adaptive crossover");
        return {};
    }

    // Validate parents first
    std::vector<std::shared_ptr<Individual>> validParents;
    for (const auto& parent : parents) {
        if (parent && representation->isValid(parent, instance)) {
            validParents.push_back(parent);
        }
    }

    if (validParents.size() < 2) {
        LOG_WARNING("Not enough valid parents for crossover");
        return parents;  // Return original parents as fallback
    }

    // Select crossover operator based on performance
    auto crossoverOp = selectCrossover();
    if (!crossoverOp) {
        LOG_WARNING("No crossover operator selected");
        return validParents;
    }

    // Perform crossover with valid parents
    auto offspring = crossoverOp->crossover(validParents, instance, representation);
    
    // If crossover produced no offspring (technical error), return valid parents
    if (offspring.empty()) {
        LOG_WARNING("Crossover produced no offspring - technical error");
        return validParents;
    }
    
    // Keep offspring regardless of validation status
    // Let fitness function handle penalties for mismatches
    return offspring;
}

void AdaptiveCrossover::updateFeedback(double currentBestFitness) {
    LOG_DEBUG("Entering updateFeedback, generation: " + std::to_string(generationCount));
    
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
    
    // Log cache and population state for debugging
    if (m_config.getCache()) {
        auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
        if (cache) {
            const auto& population = cache->getCurrentPopulation();
            LOG_DEBUG("Cache and population state - Cache: valid, Population size: " + 
                     std::to_string(population.size()) + ", Generation: " + 
                     std::to_string(generationCount));
        } else {
            LOG_DEBUG("Cache and population state - Cache: invalid cast");
        }
    } else {
        LOG_DEBUG("Cache and population state - Cache: null");
    }
}

