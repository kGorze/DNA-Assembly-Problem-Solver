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
    // Get singleton instance
    auto& config = GAConfig::getInstance();
    
    try {
        config.loadFromFile("config.cfg");
    } catch (const std::exception& e) {
        std::cerr << "Failed to load GA configuration: " << e.what() << std::endl;
        return;
    }
    
    std::ofstream csvFile(outputFile);
    
    // Write CSV header
    csvFile << "Inertia,AdaptationInterval,MinTrials,MinProb,RunNumber,"
            << "AvgFitness,BestFitness,ConvergenceGen,ExecutionTime,"
            << "OX_Usage,ERX_Usage,PMX_Usage,"
            << "OX_Success,ERX_Success,PMX_Success\n";

    // Define parameter ranges for the benchmark
    std::vector<double> inertiaValues = {0.1, 0.3, 0.5, 0.7, 0.9};
    std::vector<int> adaptationIntervals = {5, 10, 20};
    std::vector<int> minTrialsValues = {3, 5, 10};
    std::vector<double> minProbValues = {0.1, 0.2, 0.3};
    int runsPerConfig = 5;

    for (double inertia : inertiaValues) {
        for (int interval : adaptationIntervals) {
            for (int trials : minTrialsValues) {
                for (double minProb : minProbValues) {
                    for (int run = 0; run < runsPerConfig; run++) {
                        // Configure adaptive crossover with current parameters
                        auto crossover = std::make_shared<AdaptiveCrossover>();
                        crossover->setParameters(inertia, interval, trials, minProb);
                        
                        // Run GA with this configuration
                        auto startTime = std::chrono::high_resolution_clock::now();
                        
                        // Create GA components
                        auto cache = std::make_shared<SimplePopulationCache>();
                        config.setCache(cache);

                        auto ga = GeneticAlgorithm(
                            std::make_shared<PermutationRepresentation>(),
                            std::make_shared<TournamentSelection>(config, cache),
                            crossover,
                            std::make_shared<PointMutation>(config.getMutationRate()),
                            std::make_shared<PartialReplacement>(config.getReplacementRatio(), cache),
                            std::make_shared<OptimizedGraphBasedFitness>(),
                            std::make_shared<MaxGenerationsStopping>(config),
                            cache,
                            config
                        );
                        
                        ga.run(instance);
                        
                        auto endTime = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
                                      (endTime - startTime);

                        // Get metrics from the run
                        RunMetrics metrics = crossover->getMetrics();
                        metrics.executionTime = duration.count();

                        // Write results to CSV
                        csvFile << std::fixed << std::setprecision(4)
                                << inertia << ","
                                << interval << ","
                                << trials << ","
                                << minProb << ","
                                << run << ","
                                << metrics.avgFitness << ","
                                << metrics.bestFitness << ","
                                << metrics.convergenceTime << ","
                                << metrics.executionTime;

                        // Write operator usage rates
                        for (double rate : metrics.operatorUsageRates) {
                            csvFile << "," << rate;
                        }

                        // Write operator success rates
                        for (double rate : metrics.operatorSuccessRates) {
                            csvFile << "," << rate;
                        }
                        csvFile << "\n";
                    }
                }
            }
        }
    }
    
    csvFile.close();
    std::cout << "Benchmark results saved to: " << outputFile << std::endl;
}
