//
// Created by konrad_guest on 11/01/2025.
// SMART

#ifndef ADAPTIVE_CROSSOVER_BENCHMARK_H
#define ADAPTIVE_CROSSOVER_BENCHMARK_H

#include <vector>
#include <string>
#include "metaheuristics/adaptive_crossover.h"
#include "generator/dna_generator.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "metaheuristics/genetic_algorithm.h"



struct BenchmarkResults {
    double inertia;
    int adaptationInterval;
    int minTrials;
    double minProb;
    RunMetrics metrics;
};

class AdaptiveCrossoverBenchmark {
private:
    struct BenchmarkConfig {
        std::vector<double> inertiaValues = {0.5, 0.7, 0.9};
        std::vector<int> adaptationIntervals = {10, 20, 30};
        std::vector<int> minTrials = {3, 5, 7};
        std::vector<double> minProbs = {0.05, 0.1, 0.15};
        int runsPerConfig = 5;
        int maxGenerations = 1000;
    };

    std::vector<BenchmarkResults> results;
    std::string outputFile;

public:
    AdaptiveCrossoverBenchmark(const std::string& outputPath = "adaptive_crossover_benchmark.csv");

    void runBenchmark(const DNAInstance& instance);
};

#endif //ADAPTIVE_CROSSOVER_BENCHMARK_H
