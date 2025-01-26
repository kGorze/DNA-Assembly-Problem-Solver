//
// Created by konrad_guest on 08/01/2025.
// SMART


#include "../../include/benchmark/crossover_benchmark.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/metaheuristics/adaptive_crossover.h"
#include "../../include/utils/logging.h"
#include <chrono>
#include <random>
#include <algorithm>
#include <numeric>

CrossoverBenchmark::CrossoverBenchmark(double crossoverRate)
    : m_crossoverRate(crossoverRate)
{
    // Create adaptive crossover with its operators
    auto adaptiveCrossover = std::make_shared<AdaptiveCrossover>();
    m_crossovers.push_back(adaptiveCrossover);
}

void CrossoverBenchmark::runBenchmark(
    const std::vector<std::shared_ptr<Individual>>& population,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation,
    int iterations)
{
    LOG_INFO("Starting adaptive crossover benchmark with " + std::to_string(iterations) + " iterations");
    
    BenchmarkResult result;
    result.totalTime = 0;
    result.successRate = 0;
    result.avgOffspringFitness = 0;
    
    int validOffspring = 0;
    double totalFitness = 0;
    
    for (int i = 0; i < iterations; i++) {
        // Select random parents
        std::vector<std::shared_ptr<Individual>> parents;
        parents.push_back(population[rand() % population.size()]);
        parents.push_back(population[rand() % population.size()]);
        
        // Measure crossover time
        auto start = std::chrono::high_resolution_clock::now();
        auto offspring = m_crossovers[0]->crossover(parents, instance, representation);
        auto end = std::chrono::high_resolution_clock::now();
        
        result.totalTime += std::chrono::duration<double>(end - start).count();
        
        // Calculate success rate and fitness
        for (const auto& child : offspring) {
            if (child && representation->isValid(child, instance)) {
                validOffspring++;
                totalFitness += child->getFitness();
            }
        }
    }
    
    // Calculate final metrics
    result.successRate = static_cast<double>(validOffspring) / (iterations * 2);
    result.avgOffspringFitness = validOffspring > 0 ? totalFitness / validOffspring : 0;
    
    m_results.push_back(result);
    
    LOG_INFO("Adaptive crossover benchmark completed");
}

void CrossoverBenchmark::printResults() const {
    std::cout << "\nAdaptive Crossover Benchmark Results:\n";
    std::cout << "-----------------------------------\n";
    std::cout << "Average time: " << m_results[0].totalTime << " ms\n";
    std::cout << "Success rate: " << m_results[0].successRate * 100 << "%\n";
    std::cout << "Avg offspring fitness: " << m_results[0].avgOffspringFitness << "\n\n";
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
    
    // Create components with permutation representation and optimized graph-based fitness
    auto cache = std::make_shared<SimplePopulationCache>();
    config.setCache(cache);
    
    auto representation = std::make_shared<PermutationRepresentation>();
    auto selection = config.getSelection();
    auto mutation = std::make_shared<OnePointMutation>();
    auto replacement = config.getReplacement();
    auto fitness = std::make_shared<OptimizedGraphBasedFitness>();
    auto stopping = config.getStopping();
    
    // Create and run GA with adaptive crossover
    GeneticAlgorithm ga(
        representation,
        selection,
        crossover,
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
