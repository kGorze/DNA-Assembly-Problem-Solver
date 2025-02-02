//
// Created by konrad_guest on 10/01/2025. - test already in the suite.
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
#include "utils/random.h"

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
    int convergenceGeneration;
    std::vector<std::vector<double>> operatorUsageHistory;
    std::vector<std::vector<double>> operatorSuccessHistory;
    std::vector<double> fitnessHistory;
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
public:
    explicit AdaptiveCrossover(const GAConfig& config);
    ~AdaptiveCrossover() override = default;

    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;

    void updateFeedback(double fitness);
    void logMetrics() const;
    void setGeneration(int gen) { generationCount = gen; }
    void setParameters(double inertia, int adaptInterval, int minTrials, double minProb);
    RunMetrics getMetrics() const;

    // Static constants
    inline static constexpr double INERTIA = 0.5;
    inline static constexpr int ADAPTATION_INTERVAL = 10;
    inline static constexpr int MIN_TRIALS = 5;
    inline static constexpr double MIN_PROB = 0.1;
    inline static constexpr double EPSILON = 1e-6;

private:
    const GAConfig& m_config;
    std::unique_ptr<Random> m_random;
    std::vector<std::shared_ptr<ICrossover>> m_crossovers;
    std::vector<double> m_crossoverPerformance;
    std::vector<int> m_crossoverUsage;
    std::vector<double> m_crossoverProbabilities;
    std::vector<int> m_crossoverSuccesses;
    std::vector<int> m_crossoverTrials;
    int m_currentCrossoverIndex;
    double m_lastDiversityMeasure;
    double previousBestFitness;
    double bestSeenFitness;
    int generationCount;
    RunMetrics metrics;

    std::shared_ptr<ICrossover> selectCrossover();
    void updatePerformance(bool improved);
    void adjustProbabilities();
    void logDiversityMetrics(
        const std::vector<std::shared_ptr<Individual>>& population,
        std::shared_ptr<IRepresentation> representation,
        RunMetrics& metrics);
    double calculateAverageDistance(
        const std::vector<std::shared_ptr<Individual>>& population) const;
    void updateMetrics(
        const std::vector<std::shared_ptr<Individual>>& offspring,
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation);
};
