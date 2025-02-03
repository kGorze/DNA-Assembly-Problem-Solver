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
#include "../utils/logging.h"
#include "../dna/dna_instance.h"

#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <mutex>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <queue>
#include <random>
#include <chrono>

/**
 * Simple fitness function that counts the number of k-mers
 * shared with the original DNA spectrum.
 */
class SimpleFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

private:
    double calculateConnectivityScore(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const;
    double calculateSpectrumCoverageScore(const std::vector<char>& dna, const DNAInstance& instance) const;
    double calculateLengthPenalty(int actualLength, int targetLength) const;
};

/**
 * Better fitness function that uses more sophisticated metrics.
 */
class BetterFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

private:
    double calculateConnectivityScore(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const;
    double calculateSpectrumCoverageScore(const std::vector<char>& dna, const DNAInstance& instance) const;
    double calculateLengthPenalty(int actualLength, int targetLength) const;
    double calculateSequenceSimilarity(const std::vector<char>& dna, const DNAInstance& instance) const;
};

/**
 * Smith-Waterman based fitness function.
 */
class SmithWatermanFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

private:
    int smithWaterman(const std::string& seq1, const std::string& seq2) const;
};

/**
 * Optimized graph-based fitness function.
 */
class OptimizedGraphBasedFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

protected:
    // Protected methods for testing
    void buildGraphForTest(const std::vector<int>& genes, const DNAInstance& instance) const {
        buildGraph(genes, instance);
    }

    bool testNodesConnection(int node1, int node2, const DNAInstance& instance) const {
        return areNodesConnected(node1, node2, instance);
    }

    double getConnectivityScore() const {
        return calculateConnectivityScore();
    }

    double getLengthPenalty(const std::vector<int>& genes, const DNAInstance& instance) const {
        return calculateLengthPenalty(genes, instance);
    }

    double getSpectrumCoverageScore(const std::vector<int>& genes, const DNAInstance& instance) const {
        return calculateSpectrumCoverageScore(genes, instance);
    }

    const std::vector<std::vector<bool>>& getAdjacencyMatrix() const {
        return adjacencyMatrix;
    }

private:
    void ensureBufferSizes(size_t size) const;
    void buildGraph(const std::vector<int>& genes, const DNAInstance& instance) const;
    bool areNodesConnected(int node1, int node2, const DNAInstance& instance) const;
    double calculateConnectivityScore() const;
    double calculateSpectrumCoverageScore(const std::vector<int>& genes, const DNAInstance& instance) const;
    double calculateLengthPenalty(const std::vector<int>& genes, const DNAInstance& instance) const;
    void cleanupCache() const;

    mutable std::vector<std::vector<bool>> adjacencyMatrix;
    mutable std::vector<std::vector<int>> distanceMatrix;
    mutable std::mutex graphCacheMutex;
};
