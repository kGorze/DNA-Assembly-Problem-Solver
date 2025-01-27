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
#include <set>

// Define DiversityMetrics before it's used
struct DiversityMetrics {
    double avgHammingDistance;
    double uniqueSolutionsRatio;
    double uniqueDNASequencesRatio;
};

struct RunMetrics {
    double avgFitness;
    double bestFitness;
    std::vector<double> operatorUsageRates;
    std::vector<double> operatorSuccessRates;
    double convergenceTime;  // in generations
    double executionTime;    // in milliseconds
    std::vector<DiversityMetrics> diversityHistory;
};

// Define CrossoverPerformance before it's used
struct CrossoverPerformance {
    std::shared_ptr<ICrossover> crossover;
    double successRate;
    double recentSuccessRate;
    int trials;
    int successes;
    int usageCount;
    int successCount;
    int recentUsageCount;
    int recentSuccessCount;

    explicit CrossoverPerformance(std::shared_ptr<ICrossover> op) 
        : crossover(op)
        , successRate(0.0)
        , recentSuccessRate(0.0)
        , trials(0)
        , successes(0)
        , usageCount(0)
        , successCount(0)
        , recentUsageCount(0)
        , recentSuccessCount(0) {}
};

class AdaptiveCrossover : public ICrossover {
private:
    // Reorder members to match initialization order
    GAConfig m_config;
    DNAInstance m_instance;
    std::unique_ptr<Random> m_random;
    double m_lastDiversityMeasure;
    
    // Add member variables for crossover operators
    std::shared_ptr<ICrossover> m_orderCrossover;
    std::shared_ptr<ICrossover> m_edgeRecombination;
    std::shared_ptr<ICrossover> m_pmxCrossover;
    std::shared_ptr<ICrossover> m_dnaAlignmentCrossover;
    
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
        std::vector<DiversityMetrics> diversityHistory;
        int convergenceGeneration;
        double bestFitness;
        double avgFitness;
    } metrics;

    std::vector<CrossoverPerformance> crossovers;
    
    void updatePerformance(bool improved);
    void adjustProbabilities();
    std::shared_ptr<ICrossover> selectCrossover();
    
    void logDiversityMetrics();
    double calculateAverageDistance(const std::vector<std::shared_ptr<Individual>>& population) const;

public:
    AdaptiveCrossover(const GAConfig& config, const DNAInstance& instance);
    
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
        
    void setParameters(double inertia, int adaptInterval, int minTrials, double minProb);
    void setGeneration(int gen) { generationCount = gen; }
    void updateFeedback(double currentBestFitness);
    RunMetrics getMetrics() const;
};
