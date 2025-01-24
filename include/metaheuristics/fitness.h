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

#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>

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

class OptimizedGraphBasedFitness : public IFitness {
public:
    struct Edge {
        int to;
        int weight;
        Edge(int t, int w) : to(t), weight(w) {}
    };

    struct PreprocessedEdge {
        int to;
        int weight;
        bool exists;
        PreprocessedEdge() : to(-1), weight(0), exists(false) {}
        PreprocessedEdge(int t, int w, bool e) : to(t), weight(w), exists(e) {}
    };

    struct PathAnalysis {
        int edgesWeight1;
        int edgesWeight2or3;
        int uniqueNodesUsed;
        int repeatNodeUsages;
        PathAnalysis() : edgesWeight1(0), edgesWeight2or3(0), uniqueNodesUsed(0), repeatNodeUsages(0) {}
        PathAnalysis(int w1, int w23, int unique, int repeat) 
            : edgesWeight1(w1), edgesWeight2or3(w23), uniqueNodesUsed(unique), repeatNodeUsages(repeat) {}
    };

    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

private:
    mutable std::vector<int> nodeUsageBuffer;
    mutable std::unordered_map<std::string, std::vector<std::vector<Edge>>> graphCache;

    void initBuffers(size_t size) const;
    std::string createCacheKey(const std::vector<std::string>& spectrum, int k) const;
    std::vector<std::vector<Edge>> buildSpectrumGraph(const std::vector<std::string>& spectrum, int k) const;
    int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;
    PathAnalysis analyzePath(const std::vector<int>& path, const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix) const;
    const std::vector<int>& permutationToPath(std::shared_ptr<std::vector<int>> individual) const;
    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(
        const std::vector<std::vector<Edge>>& graph
    ) const;
};

class Fitness : public IFitness, public std::enable_shared_from_this<const IFitness> {
private:
    std::shared_ptr<IPopulationCache> m_cache;

protected:
    virtual double calculateDNAFitness(const std::vector<char>& dna, const DNAInstance& instance) const;
    int levenshteinDistance(const std::string& s1, const std::string& s2) const;

public:
    explicit Fitness(std::shared_ptr<IPopulationCache> cache = nullptr) : m_cache(cache) {}
    
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;
};
