//
// Created by konrad_guest on 10/01/2025.
// SMART

#pragma once

#include "../interfaces/i_crossover.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "crossover_impl.h"
#include <vector>
#include <memory>
#include <random>
#include <string>

struct RunMetrics {
    double avgFitness;
    double bestFitness;
    std::vector<double> operatorUsageRates;
    std::vector<double> operatorSuccessRates;
    double convergenceTime;  // in generations
    double executionTime;    // in milliseconds
};

class AdaptiveCrossover : public ICrossover {
private:
    struct CrossoverPerformance {
        std::shared_ptr<ICrossover> crossover;
        double successRate;
        double recentSuccessRate;  
        int usageCount;
        int successCount;
        int recentUsageCount;    
        int recentSuccessCount;   
        
        CrossoverPerformance(std::shared_ptr<ICrossover> c) 
            : crossover(c), successRate(0.33), recentSuccessRate(0.33),
              usageCount(0), successCount(0), 
              recentUsageCount(0), recentSuccessCount(0) {}
    };
    
    std::vector<CrossoverPerformance> crossovers;
    double previousBestFitness;
    double bestSeenFitness;       
    int currentCrossoverIndex;
    int generationCount;
    
    // Configurable parameters
    double INERTIA;
    int ADAPTATION_INTERVAL;
    int MIN_TRIALS;
    double MIN_PROB;
    const double EPSILON = 1e-6;   
    
    struct Metrics {
        std::vector<double> fitnessHistory;
        std::vector<std::vector<double>> operatorUsageHistory;
        std::vector<std::vector<double>> operatorSuccessHistory;
        int convergenceGeneration;
        double bestFitness;
        double avgFitness;
    } metrics;

    // Add missing member variables
    GAConfig m_config;
    std::unique_ptr<Random> m_random;
    
    // Store crossover operators to ensure proper lifetime
    std::shared_ptr<ICrossover> m_orderCrossover;
    std::shared_ptr<ICrossover> m_edgeRecombination;
    std::shared_ptr<ICrossover> m_pmxCrossover;
    std::shared_ptr<ICrossover> m_dnaAlignmentCrossover;
    
    void updatePerformance(bool improved);
    void adjustProbabilities();
    std::shared_ptr<ICrossover> selectCrossover();
    
public:
    AdaptiveCrossover();
    explicit AdaptiveCrossover(const GAConfig& config);
    
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
        
    void updateFeedback(double currentBestFitness);
    void setParameters(double inertia, int adaptInterval, int minTrials, double minProb);
    RunMetrics getMetrics() const;
};
