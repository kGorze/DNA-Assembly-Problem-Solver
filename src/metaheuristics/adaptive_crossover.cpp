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
    
    // Log diversity metrics every generation
    if (m_config.getCache()) {
        auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
        if (cache && !cache->getCurrentPopulation().empty()) {
            LOG_DEBUG("Logging diversity metrics for generation " + std::to_string(generationCount));
            logDiversityMetrics();
        }
    }
    
    auto& current = crossovers[currentCrossoverIndex];
    current.usageCount++;
    current.recentUsageCount++;
    
    // Calculate improvement magnitude
    double improvementMagnitude = std::max(0.0, previousBestFitness - bestSeenFitness);
    bool significantImprovement = improved || (improvementMagnitude > EPSILON);
    
    // Decay old success counts
    for (auto& crossover : crossovers) {
        crossover.successCount = static_cast<int>(crossover.successCount * 0.95);  // 5% decay per generation
        crossover.recentSuccessCount = static_cast<int>(crossover.recentSuccessCount * 0.9);  // 10% decay for recent
    }
    
    if (significantImprovement) {
        current.successCount++;
        current.recentSuccessCount++;
        
        // Scale success rate by improvement magnitude
        double improvementFactor = 1.0 + std::min(1.0, improvementMagnitude * 2.0);
        bestSeenFitness = previousBestFitness;
        
        // Update recent success rate with improvement scaling
        if (current.recentUsageCount > 0) {
            current.recentSuccessRate = (static_cast<double>(current.recentSuccessCount) / 
                                       current.recentUsageCount) * improvementFactor;
        }
        
        // Add diversity bonus if this operator produces more diverse offspring
        if (m_config.getCache()) {
            auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
            if (cache) {
                const auto& population = cache->getCurrentPopulation();
                double avgDistance = calculateAverageDistance(population);
                if (avgDistance > m_lastDiversityMeasure) {
                    current.recentSuccessRate *= 1.2;  // 20% bonus for diversity increase
                }
                m_lastDiversityMeasure = avgDistance;
            }
        }
    } else {
        // Penalize lack of improvement more severely if diversity is low
        double diversityPenalty = 1.0;
        if (m_config.getCache()) {
            auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
            if (cache) {
                const auto& population = cache->getCurrentPopulation();
                double avgDistance = calculateAverageDistance(population);
                if (avgDistance < 0.2) {  // Low diversity threshold
                    diversityPenalty = 0.7;  // 30% extra penalty
                }
            }
        }
        
        if (current.recentUsageCount > 0) {
            current.recentSuccessRate = (static_cast<double>(current.recentSuccessCount) / 
                                       current.recentUsageCount) * 0.9 * diversityPenalty;
        }
    }
    
    // Update overall success rate with inertia and diversity consideration
    if (current.usageCount > 0) {
        double historicalRate = static_cast<double>(current.successCount) / current.usageCount;
        current.successRate = (INERTIA * current.successRate) + 
                            ((1.0 - INERTIA) * current.recentSuccessRate);
        
        // Blend with historical rate to prevent over-specialization
        current.successRate = (0.7 * current.successRate) + (0.3 * historicalRate);
    }
    
    generationCount++;  // Move increment to end of updatePerformance
    LOG_DEBUG("Generation count after increment: " + std::to_string(generationCount));
    
    // Periodically adjust probabilities and log diversity
    if (generationCount % ADAPTATION_INTERVAL == 0) {
        LOG_DEBUG("Reached adaptation interval at generation " + std::to_string(generationCount));
        adjustProbabilities();
        // Log diversity metrics after population has been updated
        if (m_config.getCache()) {
            auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
            if (cache && !cache->getCurrentPopulation().empty()) {
                LOG_DEBUG("Logging diversity metrics at generation " + std::to_string(generationCount));
                logDiversityMetrics();
            } else {
                LOG_DEBUG("Cache is null or population is empty at generation " + std::to_string(generationCount));
            }
        } else {
            LOG_DEBUG("No cache available at generation " + std::to_string(generationCount));
        }
    } else {
        LOG_DEBUG("Not at adaptation interval: " + std::to_string(generationCount) + " % " + 
                  std::to_string(ADAPTATION_INTERVAL) + " = " + 
                  std::to_string(generationCount % ADAPTATION_INTERVAL));
    }
}

void AdaptiveCrossover::adjustProbabilities() {
    // Debug logging for adaptation interval
    LOG_DEBUG("Adjusting probabilities at generation " + std::to_string(generationCount));
    
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
    
    // Log diversity metrics after adjusting probabilities
    if (m_config.getCache()) {
        auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
        if (cache && !cache->getCurrentPopulation().empty()) {
            LOG_DEBUG("Logging diversity metrics after probability adjustment");
            logDiversityMetrics();
        } else {
            LOG_DEBUG("Cannot log diversity metrics after probability adjustment - invalid cache or empty population");
        }
    } else {
        LOG_DEBUG("Cannot log diversity metrics after probability adjustment - no cache");
    }
}

void AdaptiveCrossover::logDiversityMetrics() {
    LOG_DEBUG("Attempting to log diversity metrics at generation " + std::to_string(generationCount));
    
    if (!m_config.getCache()) {
        LOG_DEBUG("No cache available for diversity metrics");
        return;
    }
    
    auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache());
    if (!cache) {
        LOG_DEBUG("Cache cast failed for diversity metrics");
        return;
    }
    
    const auto& population = cache->getCurrentPopulation();
    if (population.empty()) {
        LOG_DEBUG("Population is empty, cannot calculate diversity metrics");
        return;
    }

    // Calculate average Hamming distance
    double avgDistance = calculateAverageDistance(population);
    
    // Calculate unique solutions ratio
    std::set<std::vector<int>> uniqueGenes;
    for (const auto& individual : population) {
        if (individual && !individual->getGenes().empty()) {
            uniqueGenes.insert(individual->getGenes());
        }
    }
    double uniqueSolutionsRatio = static_cast<double>(uniqueGenes.size()) / population.size();
    
    // Calculate unique DNA sequences ratio with safety checks
    std::set<std::string> uniqueDNA;
    auto representation = m_config.getRepresentation();
    if (representation) {
        for (const auto& individual : population) {
            if (!individual || individual->getGenes().empty()) {
                continue;
            }
            
            // Validate genes before conversion
            bool validGenes = true;
            const auto& genes = individual->getGenes();
            const auto& spectrum = m_instance.getSpectrum();
            
            for (int gene : genes) {
                if (gene < 0 || static_cast<size_t>(gene) >= spectrum.size()) {
                    validGenes = false;
                    break;
                }
            }
            
            if (!validGenes) {
                LOG_DEBUG("Skipping individual with invalid genes");
                continue;
            }
            
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
    } else {
        LOG_DEBUG("No representation available for DNA conversion");
    }
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

