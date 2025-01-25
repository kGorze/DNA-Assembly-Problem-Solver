//
// Created by konrad_guest on 07/01/2025.
// SMART

#pragma once

#include "../interfaces/i_representation.h"
#include "../interfaces/i_fitness.h"
#include "../interfaces/i_population_cache.h"
#include "../generator/dna_generator.h"
#include "../utils/utility_functions.h"
#include "../metaheuristics/path_analyzer.h"
#include "../metaheuristics/preprocessed_edge.h"
#include "../metaheuristics/individual.h"
#include "../utils/logger.h"
#include "../dna/dna_instance.h"

#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <mutex>
#include <cmath>
#include <limits>

/**
 * Bardzo prosta funkcja fitness – liczy liczbę k-merów
 * wspólnych z oryginalnym spektrum DNA.
 */
class SimpleFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;
};

/**
 * Ulepszona wersja fitness – oprócz liczby trafionych k-merów
 * uwzględnia odległość Levenshteina do oryginalnego DNA.
 */
class BetterFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;
};

/**
 * Przykład użycia algorytmu Smith-Waterman do oceny podobieństwa.
 */
class SmithWatermanFitness : public IFitness {
public:
    static constexpr int MATCH_SCORE = 2;
    static constexpr int MISMATCH_SCORE = -1;
    static constexpr int GAP_PENALTY = -2;

    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

private:
    int smithWaterman(const std::string& seq1, const std::string& seq2) const;
};

class IFitness {
public:
    virtual ~IFitness() = default;
    virtual double evaluate(const std::shared_ptr<Individual>& individual,
                          const DNAInstance& instance) = 0;
};

class Fitness : public IFitness {
protected:
    std::shared_ptr<IPopulationCache> m_cache;
    
    virtual double calculateDNAFitness(const std::string& dna, const DNAInstance& instance) const = 0;
    
    double levenshteinDistance(const std::string& s1, const std::string& s2) const {
        if (s1.empty()) return s2.length();
        if (s2.empty()) return s1.length();
        
        std::vector<std::vector<double>> dp(s1.length() + 1, 
                                          std::vector<double>(s2.length() + 1));
        
        for (size_t i = 0; i <= s1.length(); i++) dp[i][0] = i;
        for (size_t j = 0; j <= s2.length(); j++) dp[0][j] = j;
        
        for (size_t i = 1; i <= s1.length(); i++) {
            for (size_t j = 1; j <= s2.length(); j++) {
                dp[i][j] = std::min({
                    dp[i-1][j] + 1,
                    dp[i][j-1] + 1,
                    dp[i-1][j-1] + (s1[i-1] != s2[j-1])
                });
            }
        }
        
        return dp[s1.length()][s2.length()];
    }
    
public:
    explicit Fitness(std::shared_ptr<IPopulationCache> cache) : m_cache(std::move(cache)) {}
    virtual ~Fitness() = default;
    
    double evaluate(const std::shared_ptr<Individual>& individual,
                   const DNAInstance& instance) override {
        if (!individual) {
            LOG_ERROR("Cannot evaluate null individual");
            return -std::numeric_limits<double>::infinity();
        }
        
        if (!m_cache) {
            LOG_ERROR("Cache not initialized");
            return -std::numeric_limits<double>::infinity();
        }
        
        try {
            return m_cache->getOrCalculateFitness(individual, instance);
        } catch (const std::exception& e) {
            LOG_ERROR("Fitness evaluation failed: " + std::string(e.what()));
            return -std::numeric_limits<double>::infinity();
        }
    }
};

struct PreprocessedEdge {
    bool exists = false;
    double weight = 0.0;
    
    PreprocessedEdge() = default;
    PreprocessedEdge(bool e, double w) : exists(e), weight(w) {}
};

class OptimizedGraphBasedFitness : public IFitness {
private:
    static constexpr size_t maxCacheSize = 10000;
    std::unordered_map<std::string, double> cache;
    std::mutex graphCacheMutex;
    
    // Pre-allocated buffers for performance
    std::vector<int> nodeUsageBuffer;
    std::vector<bool> coveredBuffer;
    std::vector<std::vector<PreprocessedEdge>> adjacencyMatrixBuffer;
    
    void ensureBufferSizes(size_t size);
    void buildGraph(const std::vector<int>& genes, const DNAInstance& instance);
    bool areNodesConnected(int node1, int node2, const DNAInstance& instance) const;
    double calculateConnectivityScore() const;
    double calculateSpectrumCoverageScore(const std::vector<int>& genes, const DNAInstance& instance) const;
    double calculateLengthPenalty(const std::vector<int>& genes, const DNAInstance& instance) const;
    void cleanupCache();
    
public:
    OptimizedGraphBasedFitness() = default;
    ~OptimizedGraphBasedFitness() override = default;
    
    double evaluate(const std::shared_ptr<Individual>& individual,
                   const DNAInstance& instance) override;
};
