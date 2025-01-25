//
// Created by konrad_guest on 08/01/2025.
// SMART


#include "../../include/benchmark/crossover_benchmark.h"
#include "../../include/metaheuristics/representation_impl.h"
#include "../../include/metaheuristics/selection_impl.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/replacement_impl.h"
#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
#include "../../include/metaheuristics/population_cache_impl.h"
#include "utils/logging.h"

void CrossoverBenchmark::runBenchmark(const DNAInstance& instance) {
    // Create config for crossover rate
    GAConfig config;
    if (!config.loadFromFile("config.cfg")) {
        throw std::runtime_error("Failed to load configuration");
    }
    double crossoverRate = config.getCrossoverProbability();
    
    // Create crossover operators with proper rate
    auto onePoint = std::make_shared<OnePointCrossover>(crossoverRate);
    auto order = std::make_shared<OrderCrossover>();
    auto edge = std::make_shared<EdgeRecombination>();
    auto pmx = std::make_shared<PMXCrossover>();
    auto distPreserving = std::make_shared<DistancePreservingCrossover>();
    
    // Create crossover operators to benchmark
    std::vector<std::pair<std::string, std::shared_ptr<ICrossover>>> crossovers = {
        {"OnePoint", onePoint},
        {"OrderCrossover", order},
        {"EdgeRecombination", edge},
        {"PMXCrossover", pmx},
        {"DistancePreserving", distPreserving}
    };
    
    // Create representation
    auto representation = std::make_shared<PermutationRepresentation>();
    
    // Run benchmarks
    for (const auto& [name, crossover] : crossovers) {
        LOG_INFO("Running benchmark for " + name);
        
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
        for (int i = 0; i < numTrials; i++) {
            auto offspring = crossover->crossover(parents, instance, representation);
            validOffspring += offspring.size();
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        
        // Log results
        LOG_INFO(name + " results:");
        LOG_INFO("  Average time per crossover: " + 
                 std::to_string(duration.count() / (double)numTrials) + " microseconds");
        LOG_INFO("  Average offspring per crossover: " + 
                 std::to_string(validOffspring / (double)numTrials));
    }
}

int CrossoverBenchmark::runOneGA(
    const DNAInstance& instance,
    std::shared_ptr<ICrossover> crossover,
    const std::string& name)
{
    // Create local config for this run
    GAConfig config;
    if (!config.loadFromFile("config.cfg")) {
        throw std::runtime_error("Failed to load configuration");
    }
    
    // Create components
    auto cache = std::make_shared<SimplePopulationCache>();
    config.setCache(cache);
    
    auto representation = config.getRepresentation();
    auto selection = config.getSelection();
    auto mutation = config.getMutation();
    auto replacement = config.getReplacement();
    auto fitness = config.getFitness();
    auto stopping = config.getStopping();
    
    // Create and run GA
    GeneticAlgorithm ga(
        representation,
        selection,
        crossover,  // Use the provided crossover operator
        mutation,
        replacement,
        fitness,
        stopping,
        cache,
        config
    );
    
    ga.run(instance);
    return ga.getBestFitness();
}
