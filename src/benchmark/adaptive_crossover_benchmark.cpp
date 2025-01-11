//
// Created by konrad_guest on 11/01/2025.
// SMART

#include "benchmark/adaptive_crossover_benchmark.h"
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
    BenchmarkConfig config;
    std::ofstream csvFile(outputFile);
    
    // Write CSV header
    csvFile << "Inertia,AdaptationInterval,MinTrials,MinProb,RunNumber,"
            << "AvgFitness,BestFitness,ConvergenceGen,ExecutionTime,"
            << "OX_Usage,ERX_Usage,PMX_Usage,"
            << "OX_Success,ERX_Success,PMX_Success\n";

    for (double inertia : config.inertiaValues) {
        for (int interval : config.adaptationIntervals) {
            for (int trials : config.minTrials) {
                for (double minProb : config.minProbs) {
                    for (int run = 0; run < config.runsPerConfig; run++) {
                        // Configure adaptive crossover with current parameters
                        auto crossover = std::make_shared<AdaptiveCrossover>();
                        crossover->setParameters(inertia, interval, trials, minProb);
                        
                        // Run GA with this configuration
                        auto startTime = std::chrono::high_resolution_clock::now();
                        
                        // Create and run GA
                        GeneticAlgorithm ga(
                            GAConfig::getInstance().getRepresentation(),
                            GAConfig::getInstance().getSelection(),
                            crossover,
                            GAConfig::getInstance().getMutation(),
                            GAConfig::getInstance().getReplacement(),
                            std::make_shared<OptimizedGraphBasedFitness>(),
                            std::make_shared<MaxGenerationsStopping>(config.maxGenerations),
                            std::make_shared<CachedPopulation>()
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
