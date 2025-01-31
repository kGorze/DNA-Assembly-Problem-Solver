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

AdaptiveCrossover::AdaptiveCrossover(const GAConfig& config)
    : m_config(config)
    , m_random(std::make_unique<Random>())
    , m_crossovers()
    , m_crossoverPerformance()
    , m_crossoverUsage()
    , m_crossoverProbabilities()
    , m_crossoverSuccesses()
    , m_crossoverTrials()
    , m_currentCrossoverIndex(0)
    , m_lastDiversityMeasure(0.0)
    , previousBestFitness(-std::numeric_limits<double>::infinity())
    , bestSeenFitness(-std::numeric_limits<double>::infinity())
    , generationCount(0) {
    
    // Initialize crossover operators
    m_crossovers.push_back(std::make_shared<OrderCrossover>());
    m_crossovers.push_back(std::make_shared<CycleCrossover>());
    m_crossovers.push_back(std::make_shared<PMXCrossover>());
    m_crossovers.push_back(std::make_shared<EdgeRecombinationCrossover>());
    
    // Initialize performance metrics
    m_crossoverPerformance.resize(m_crossovers.size(), 0.0);
    m_crossoverUsage.resize(m_crossovers.size(), 0);
    
    // Set initial probabilities
    double initialProb = 1.0 / m_crossovers.size();
    m_crossoverProbabilities.resize(m_crossovers.size(), initialProb);
    
    // Initialize metrics
    metrics.convergenceGeneration = -1;
    metrics.bestFitness = -std::numeric_limits<double>::infinity();
    metrics.avgFitness = 0.0;
    
    LOG_DEBUG("AdaptiveCrossover initialized with ADAPTATION_INTERVAL=" + 
              std::to_string(ADAPTATION_INTERVAL));
}

void AdaptiveCrossover::setParameters(double inertia, int adaptInterval, int minTrials, double minProb) {
    // This method is now empty as the parameters are defined as constexpr
}

RunMetrics AdaptiveCrossover::getMetrics() const {
    RunMetrics result;
    result.avgFitness = metrics.avgFitness;
    result.bestFitness = metrics.bestFitness;
    result.convergenceTime = metrics.convergenceGeneration;
    
    result.operatorUsageRates.resize(m_crossovers.size());
    result.operatorSuccessRates.resize(m_crossovers.size());
    
    for (size_t i = 0; i < m_crossovers.size(); i++) {
        result.operatorUsageRates[i] = (generationCount == 0)
            ? 0.0
            : static_cast<double>(m_crossoverUsage[i]) / generationCount;

        if (m_crossoverUsage[i] > 0) {
            result.operatorSuccessRates[i] = static_cast<double>(m_crossoverPerformance[i])
                                           / m_crossoverUsage[i];
        } else {
            result.operatorSuccessRates[i] = 0.0;
        }
    }
    
    return result;
}

void AdaptiveCrossover::updatePerformance(bool improved) {
    LOG_DEBUG("Entering updatePerformance, generation: " + std::to_string(generationCount));
    
    // Update performance metrics for current crossover operator
    if (m_currentCrossoverIndex >= 0 && m_currentCrossoverIndex < static_cast<int>(m_crossovers.size())) {
        auto& currentOp = m_crossovers[m_currentCrossoverIndex];
        m_crossoverUsage[m_currentCrossoverIndex]++;
        if (improved) {
            m_crossoverPerformance[m_currentCrossoverIndex]++;
        }
        
        // Update success rates
        if (m_crossoverUsage[m_currentCrossoverIndex] > 0) {
            m_crossoverProbabilities[m_currentCrossoverIndex] = static_cast<double>(m_crossoverPerformance[m_currentCrossoverIndex])
                                                           / m_crossoverUsage[m_currentCrossoverIndex];
        }
        
        // Log operator performance
        LOG_DEBUG("Operator " + std::to_string(m_currentCrossoverIndex) + 
                 " performance - Success rate: " + std::to_string(m_crossoverProbabilities[m_currentCrossoverIndex]));
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
    for (const auto& op : m_crossovers) {
        usageRates.push_back(static_cast<double>(m_crossoverUsage[&op - &m_crossovers[0]]) / std::max(1, ADAPTATION_INTERVAL));
        successRates.push_back(m_crossoverProbabilities[&op - &m_crossovers[0]]);
    }
    metrics.operatorUsageHistory.push_back(usageRates);
    metrics.operatorSuccessHistory.push_back(successRates);
    
    // Reset recent counters if adaptation interval is reached
    if (generationCount % ADAPTATION_INTERVAL == 0) {
        for (auto& op : m_crossovers) {
            m_crossoverUsage[&op - &m_crossovers[0]] = 0;
            m_crossoverPerformance[&op - &m_crossovers[0]] = 0;
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
    
    // Calculate total success rate and find best performing operator
    double totalSuccessRate = 0.0;
    double bestSuccessRate = 0.0;
    size_t bestOperatorIdx = 0;
    
    for (size_t i = 0; i < m_crossovers.size(); i++) {
        double successRate = m_crossoverUsage[i] > 0 
            ? static_cast<double>(m_crossoverPerformance[i]) / m_crossoverUsage[i]
            : 0.0;
        totalSuccessRate += successRate;
        
        if (successRate > bestSuccessRate) {
            bestSuccessRate = successRate;
            bestOperatorIdx = i;
        }
    }
    
    // Adjust probabilities based on success rates
    if (totalSuccessRate > EPSILON) {
        // Reward successful operators more aggressively
        for (size_t i = 0; i < m_crossovers.size(); i++) {
            double successRate = m_crossoverUsage[i] > 0 
                ? static_cast<double>(m_crossoverPerformance[i]) / m_crossoverUsage[i]
                : 0.0;
            
            // Calculate new probability with stronger bias towards successful operators
            double newProb;
            if (i == bestOperatorIdx) {
                // Give best operator a higher minimum probability
                newProb = std::max(0.4, successRate / totalSuccessRate);
            } else {
                // Other operators share remaining probability based on their success
                newProb = MIN_PROB + (0.6 - MIN_PROB * (m_crossovers.size() - 1)) * 
                         (successRate / totalSuccessRate);
            }
            
            // Apply inertia to smooth probability changes
            m_crossoverProbabilities[i] = INERTIA * m_crossoverProbabilities[i] + 
                                        (1.0 - INERTIA) * newProb;
        }
    } else {
        // If no success, slightly increase probabilities of less used operators
        double totalUsage = 0.0;
        for (size_t i = 0; i < m_crossovers.size(); i++) {
            totalUsage += m_crossoverUsage[i];
        }
        
        if (totalUsage > 0) {
            for (size_t i = 0; i < m_crossovers.size(); i++) {
                double usageRatio = m_crossoverUsage[i] / totalUsage;
                // Inverse usage ratio to favor less used operators
                double newProb = MIN_PROB + (1.0 - MIN_PROB * m_crossovers.size()) * 
                               (1.0 - usageRatio);
                m_crossoverProbabilities[i] = INERTIA * m_crossoverProbabilities[i] + 
                                            (1.0 - INERTIA) * newProb;
            }
        }
    }
    
    // Normalize probabilities
    double sum = 0.0;
    for (size_t i = 0; i < m_crossovers.size(); i++) {
        sum += m_crossoverProbabilities[i];
    }
    
    if (sum > EPSILON) {
        for (size_t i = 0; i < m_crossovers.size(); i++) {
            m_crossoverProbabilities[i] /= sum;
        }
    } else {
        // Fallback to slightly biased distribution if normalization fails
        m_crossoverProbabilities[bestOperatorIdx] = 0.4;
        double remainingProb = 0.6 / (m_crossovers.size() - 1);
        for (size_t i = 0; i < m_crossovers.size(); i++) {
            if (i != bestOperatorIdx) {
                m_crossoverProbabilities[i] = remainingProb;
            }
        }
    }
    
    // Log adjusted probabilities
    std::stringstream ss;
    ss << "Adjusted probabilities:";
    for (size_t i = 0; i < m_crossovers.size(); ++i) {
        ss << "\n  Operator " << i << ": " << m_crossoverProbabilities[i];
    }
    LOG_DEBUG(ss.str());
}

void AdaptiveCrossover::logDiversityMetrics(
    const std::vector<std::shared_ptr<Individual>>& population,
    std::shared_ptr<IRepresentation> representation,
    RunMetrics& metrics) {
    
    if (population.empty()) {
        LOG_WARNING("Empty population, skipping diversity metrics");
        return;
    }
    
    // Calculate average Hamming distance with tolerance for small differences
    double totalDistance = 0.0;
    int comparisons = 0;
    
    for (size_t i = 0; i < population.size(); ++i) {
        for (size_t j = i + 1; j < population.size(); ++j) {
            if (!population[i] || !population[j]) continue;
            
            const auto& genes1 = population[i]->getGenes();
            const auto& genes2 = population[j]->getGenes();
            
            if (genes1.size() != genes2.size()) continue;
            
            int differences = 0;
            for (size_t k = 0; k < genes1.size(); ++k) {
                if (genes1[k] != genes2[k]) {  // Exact comparison, no tolerance
                    differences++;
                }
            }
            
            double distance = static_cast<double>(differences) / genes1.size();
            totalDistance += distance;
            comparisons++;
        }
    }
    
    double avgDistance = comparisons > 0 ? totalDistance / comparisons : 1.0;
    
    // Calculate unique solutions using clustering
    std::vector<std::vector<int>> clusters;  // Each cluster contains indices of similar individuals
    std::vector<bool> assigned(population.size(), false);
    
    for (size_t i = 0; i < population.size(); ++i) {
        if (!population[i] || assigned[i]) continue;
        
        // Start new cluster
        std::vector<int> cluster;
        cluster.push_back(i);
        assigned[i] = true;
        
        // Find similar individuals
        for (size_t j = i + 1; j < population.size(); ++j) {
            if (!population[j] || assigned[j]) continue;
            
            const auto& genes1 = population[i]->getGenes();
            const auto& genes2 = population[j]->getGenes();
            
            if (genes1.size() != genes2.size()) continue;
            
            // Calculate similarity
            int differences = 0;
            int matchingPositions = 0;
            for (size_t k = 0; k < genes1.size(); ++k) {
                // Count exact matches and adjacent positions
                if (genes1[k] == genes2[k]) {
                    matchingPositions++;
                } else if (k > 0 && genes1[k] == genes2[k-1] || 
                          k < genes1.size()-1 && genes1[k] == genes2[k+1]) {
                    matchingPositions++; // Count adjacent matches with half weight
                } else {
                    differences++;
                }
            }
            
            // More lenient similarity threshold - if 40% of positions match (including adjacency)
            if (static_cast<double>(matchingPositions) / genes1.size() > 0.4) {
                cluster.push_back(j);
                assigned[j] = true;
            }
        }
        clusters.push_back(cluster);
    }
    
    double uniqueSolutionsRatio = population.empty() ? 0.0 : 
                                 static_cast<double>(clusters.size()) / population.size();
    
    // Calculate unique DNA sequences using similar clustering but with different threshold
    std::vector<std::vector<int>> dnaClusters;
    std::fill(assigned.begin(), assigned.end(), false);
    
    for (size_t i = 0; i < population.size(); ++i) {
        if (!population[i] || assigned[i] || 
            !representation->isValid(population[i], *m_config.getInstance())) continue;
        
        std::vector<int> cluster;
        cluster.push_back(i);
        assigned[i] = true;
        
        std::string dna1 = representation->toDNA(population[i], *m_config.getInstance());
        if (dna1.empty()) continue;
        
        for (size_t j = i + 1; j < population.size(); ++j) {
            if (!population[j] || assigned[j] || 
                !representation->isValid(population[j], *m_config.getInstance())) continue;
                
            std::string dna2 = representation->toDNA(population[j], *m_config.getInstance());
            if (dna2.empty() || dna1.length() != dna2.length()) continue;
            
            // Calculate DNA similarity with local region matching
            int matchingRegions = 0;
            const int regionSize = 5; // Look at regions of 5 bases
            
            for (size_t k = 0; k + regionSize <= dna1.length(); k += regionSize) {
                bool regionMatches = false;
                // Look for this region anywhere within ±10 positions in the other sequence
                for (int offset = -10; offset <= 10; ++offset) {
                    if (k + offset < 0 || k + offset + regionSize > dna2.length()) continue;
                    
                    bool allMatch = true;
                    for (int r = 0; r < regionSize; ++r) {
                        if (dna1[k + r] != dna2[k + offset + r]) {
                            allMatch = false;
                            break;
                        }
                    }
                    if (allMatch) {
                        regionMatches = true;
                        break;
                    }
                }
                if (regionMatches) matchingRegions++;
            }
            
            // If 30% of regions match (allowing for shifted positions)
            double regionSimilarity = static_cast<double>(matchingRegions) / 
                                    ((dna1.length() / regionSize) + 1);
            if (regionSimilarity > 0.3) {
                cluster.push_back(j);
                assigned[j] = true;
            }
        }
        dnaClusters.push_back(cluster);
    }
    
    double uniqueDNASequencesRatio = population.empty() ? 0.0 : 
                                    static_cast<double>(dnaClusters.size()) / population.size();
    
    // Log metrics
    std::stringstream ss;
    ss << std::fixed << std::setprecision(4);
    ss << "Generation " << generationCount << " Diversity Metrics:";
    ss << "\n  Average Hamming Distance: " << avgDistance;
    ss << "\n  Unique Solutions Ratio: " << uniqueSolutionsRatio << " (" << clusters.size() << "/" << population.size() << ")";
    ss << "\n  Unique DNA Sequences Ratio: " << uniqueDNASequencesRatio << " (" << dnaClusters.size() << "/" << population.size() << ")";
    ss << "\n  Average Cluster Size: " << (population.size() / (clusters.empty() ? 1 : clusters.size()));
    
    LOG_INFO(ss.str());
    
    // Store metrics
    m_lastDiversityMeasure = avgDistance;
    DiversityMetrics diversityMetrics = {avgDistance, uniqueSolutionsRatio, uniqueDNASequencesRatio};
    metrics.diversityHistory.push_back(diversityMetrics);
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
    
    for (size_t i = 0; i < m_crossovers.size(); i++) {
        sum += m_crossoverProbabilities[i];
        if (randVal <= sum) {
            m_currentCrossoverIndex = static_cast<int>(i);
            return m_crossovers[i];
        }
    }
    
    m_currentCrossoverIndex = static_cast<int>(m_crossovers.size() - 1);
    return m_crossovers.back();
}

std::vector<std::shared_ptr<Individual>> AdaptiveCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !representation) {
        LOG_ERROR("Invalid parents or representation in adaptive crossover");
        return {};
    }
    
    // Get the shared instance from config
    auto sharedInstance = m_config.getInstance();
    if (!sharedInstance) {
        LOG_ERROR("No shared instance available in config");
        return {};
    }
    
    // Validate parents first
    std::vector<std::shared_ptr<Individual>> validParents;
    for (const auto& parent : parents) {
        if (parent && representation->isValid(parent, *sharedInstance)) {
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
        LOG_ERROR("Failed to select crossover operator");
        return {};
    }
    
    // Perform crossover
    auto offspring = crossoverOp->crossover(validParents, *sharedInstance, representation);
    
    // Update metrics
    updateMetrics(offspring, validParents, *sharedInstance, representation);
    
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

void AdaptiveCrossover::updateMetrics(
    const std::vector<std::shared_ptr<Individual>>& offspring,
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    // Implementation of updateMetrics method
    // This method should update the metrics based on the new offspring and parents
    // It should also handle the update of crossover usage and performance
}

