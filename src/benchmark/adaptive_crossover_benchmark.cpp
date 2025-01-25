//
// Created by konrad_guest on 11/01/2025.
// SMART

#include "../../include/benchmark/adaptive_crossover_benchmark.h"
#include "../../include/metaheuristics/representation_impl.h"
#include "../../include/metaheuristics/selection_impl.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/replacement_impl.h"
#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
#include "../../include/metaheuristics/population_cache_impl.h"
#include "../../include/utils/logging.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <iostream>

// Include other necessary headers for GeneticAlgorithm, GAConfig, etc.
// Make sure to include headers for GeneticAlgorithm, OptimizedGraphBasedFitness, 
// MaxGenerationsStopping, CachedPopulation, and any other dependencies.

// Constructor implementation
AdaptiveCrossoverBenchmark::AdaptiveCrossoverBenchmark(const std::string& outputPath)
    : outputFile(outputPath) {}

// runBenchmark implementation
void AdaptiveCrossoverBenchmark::runBenchmark(const DNAInstance& instance) {
    // Create and load config
    GAConfig config;
    if (!config.loadFromFile("config.cfg")) {
        throw std::runtime_error("Failed to load configuration");
    }
    
    // Create adaptive crossover with config
    auto adaptiveCrossover = std::make_shared<AdaptiveCrossover>(config);
    
    // Create representation
    auto representation = std::make_shared<PermutationRepresentation>();
    
    // Run benchmark
    LOG_INFO("Running adaptive crossover benchmark");
    
    // Create test parents
    std::vector<std::shared_ptr<std::vector<int>>> parents;
    for (int i = 0; i < 2; i++) {
        auto parent = std::make_shared<std::vector<int>>();
        *parent = representation->generateRandomSolution(instance);
        parents.push_back(parent);
    }
    
    // Measure time
    auto start = std::chrono::high_resolution_clock::now();
    
    // Run crossover multiple times
    int numTrials = 100;
    int validOffspring = 0;
    double bestFitness = -std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < numTrials; i++) {
        auto offspring = adaptiveCrossover->crossover(parents, instance, representation);
        validOffspring += offspring.size();
        
        // Simulate fitness improvement
        double currentFitness = i / (double)numTrials;
        adaptiveCrossover->updateFeedback(currentFitness);
        if (currentFitness > bestFitness) {
            bestFitness = currentFitness;
        }
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    // Get and log metrics
    auto metrics = adaptiveCrossover->getMetrics();
    
    LOG_INFO("Adaptive Crossover results:");
    LOG_INFO("  Average time per crossover: " + 
             std::to_string(duration.count() / (double)numTrials) + " microseconds");
    LOG_INFO("  Average offspring per crossover: " + 
             std::to_string(validOffspring / (double)numTrials));
    LOG_INFO("  Best fitness achieved: " + std::to_string(bestFitness));
    LOG_INFO("  Average fitness: " + std::to_string(metrics.avgFitness));
    
    // Log operator usage rates
    for (size_t i = 0; i < metrics.operatorUsageRates.size(); i++) {
        LOG_INFO("  Operator " + std::to_string(i) + " usage rate: " + 
                 std::to_string(metrics.operatorUsageRates[i]));
        LOG_INFO("  Operator " + std::to_string(i) + " success rate: " + 
                 std::to_string(metrics.operatorSuccessRates[i]));
    }
}
