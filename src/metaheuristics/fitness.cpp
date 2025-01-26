//
// Created by konrad_guest on 07/01/2025.
// SMART

#include <unordered_set>
#include <iostream>
#include "../../include/metaheuristics/fitness.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/path_analyzer.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <cmath>
#include <queue>

/* ================== SimpleFitness ================== */
double SimpleFitness::calculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!individual || !representation) {
        LOG_ERROR("Null individual or representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    try {
        std::string dnaStr = representation->toDNA(individual, instance);
        if (dnaStr.empty()) {
            LOG_ERROR("Empty DNA sequence generated from individual");
            return -std::numeric_limits<double>::infinity();
        }

        std::vector<char> dna(dnaStr.begin(), dnaStr.end());
        double connectivityScore = calculateConnectivityScore(individual, instance);
        double spectrumCoverageScore = calculateSpectrumCoverageScore(dna, instance);
        double lengthPenalty = calculateLengthPenalty(static_cast<int>(dna.size()), instance.getN());

        constexpr double CONNECTIVITY_WEIGHT = 0.4;
        constexpr double COVERAGE_WEIGHT = 0.4;
        constexpr double LENGTH_WEIGHT = 0.2;

        return (CONNECTIVITY_WEIGHT * connectivityScore) +
               (COVERAGE_WEIGHT * spectrumCoverageScore) +
               (LENGTH_WEIGHT * lengthPenalty);
    } catch (const std::exception& e) {
        LOG_ERROR("Error calculating simple fitness: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

double SimpleFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const
{
    const auto& genes = individual->getGenes();
    if (genes.empty()) return 0.0;

    double score = 0.0;
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();

    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i + 1] < 0 || genes[i + 1] >= static_cast<int>(spectrum.size())) {
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i + 1]];

        if (static_cast<int>(current.length()) >= k - 1 && 
            static_cast<int>(next.length()) >= k - 1) {
            std::string suffix = current.substr(current.length() - (k - 1));
            std::string prefix = next.substr(0, k - 1);
            if (suffix == prefix) {
                score += 1.0;
            }
        }
    }

    return score;
}

double SimpleFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const
{
    if (dna.empty()) return 0.0;

    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    
    if (static_cast<int>(dna.size()) < k) return 0.0;

    std::vector<bool> covered(spectrum.size(), false);
    double score = 0.0;

    for (size_t i = 0; i <= dna.size() - k; ++i) {
        std::string kmer(dna.begin() + i, dna.begin() + i + k);
        auto it = std::find(spectrum.begin(), spectrum.end(), kmer);
        if (it != spectrum.end()) {
            size_t index = std::distance(spectrum.begin(), it);
            if (!covered[index]) {
                covered[index] = true;
                score += 1.0;
            }
        }
    }

    return score;
}

double SimpleFitness::calculateLengthPenalty(int actualLength, int targetLength) const
{
    if (targetLength <= 0) return 0.0;
    
    double lengthDiff = std::abs(actualLength - targetLength);
    return std::exp(-lengthDiff / targetLength);
}

/* ================== BetterFitness ================== */
double BetterFitness::calculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!individual || !representation) {
        LOG_ERROR("Null individual or representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    try {
        std::string dnaStr = representation->toDNA(individual, instance);
        if (dnaStr.empty()) {
            LOG_ERROR("Empty DNA sequence generated from individual");
            return -std::numeric_limits<double>::infinity();
        }

        std::vector<char> dna(dnaStr.begin(), dnaStr.end());
        double connectivityScore = calculateConnectivityScore(individual, instance);
        double spectrumCoverageScore = calculateSpectrumCoverageScore(dna, instance);
        double lengthPenalty = calculateLengthPenalty(static_cast<int>(dna.size()), instance.getN());
        double sequenceSimilarity = calculateSequenceSimilarity(dna, instance);

        constexpr double CONNECTIVITY_WEIGHT = 0.35;
        constexpr double COVERAGE_WEIGHT = 0.35;
        constexpr double LENGTH_WEIGHT = 0.15;
        constexpr double SIMILARITY_WEIGHT = 0.15;

        return (CONNECTIVITY_WEIGHT * connectivityScore) +
               (COVERAGE_WEIGHT * spectrumCoverageScore) +
               (LENGTH_WEIGHT * lengthPenalty) +
               (SIMILARITY_WEIGHT * sequenceSimilarity);
    } catch (const std::exception& e) {
        LOG_ERROR("Error calculating better fitness: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

double BetterFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const
{
    const auto& genes = individual->getGenes();
    if (genes.empty()) return 0.0;

    double score = 0.0;
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();

    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i + 1] < 0 || genes[i + 1] >= static_cast<int>(spectrum.size())) {
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i + 1]];

        if (static_cast<int>(current.length()) >= k - 1 && 
            static_cast<int>(next.length()) >= k - 1) {
            std::string suffix = current.substr(current.length() - (k - 1));
            std::string prefix = next.substr(0, k - 1);
            if (suffix == prefix) {
                score += 1.0;
            }
        }
    }

    return score;
}

double BetterFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const
{
    if (dna.empty()) return 0.0;

    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    
    if (static_cast<int>(dna.size()) < k) return 0.0;

    std::vector<bool> covered(spectrum.size(), false);
    double score = 0.0;

    for (size_t i = 0; i <= dna.size() - k; ++i) {
        std::string kmer(dna.begin() + i, dna.begin() + i + k);
        auto it = std::find(spectrum.begin(), spectrum.end(), kmer);
        if (it != spectrum.end()) {
            size_t index = std::distance(spectrum.begin(), it);
            if (!covered[index]) {
                covered[index] = true;
                score += 1.0;
            }
        }
    }

    return score;
}

double BetterFitness::calculateLengthPenalty(int actualLength, int targetLength) const
{
    if (targetLength <= 0) return 0.0;
    
    double lengthDiff = std::abs(actualLength - targetLength);
    return std::exp(-lengthDiff / targetLength);
}

double BetterFitness::calculateSequenceSimilarity(
    const std::vector<char>& dna,
    const DNAInstance& instance) const
{
    const std::string& originalDNA = instance.getOriginalDNA();
    if (originalDNA.empty() || dna.empty()) return 0.0;

    std::string dnaStr(dna.begin(), dna.end());

    std::vector<std::vector<int>> dp(dnaStr.length() + 1, 
                                   std::vector<int>(originalDNA.length() + 1));

    for (size_t i = 0; i <= dnaStr.length(); i++) {
        dp[i][0] = i;
    }
    for (size_t j = 0; j <= originalDNA.length(); j++) {
        dp[0][j] = j;
    }

    for (size_t i = 1; i <= dnaStr.length(); i++) {
        for (size_t j = 1; j <= originalDNA.length(); j++) {
            dp[i][j] = std::min({
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,
                dp[i-1][j-1] + (dnaStr[i-1] != originalDNA[j-1])
            });
        }
    }

    double distance = dp[dnaStr.length()][originalDNA.length()];
    return 1.0 / (1.0 + distance);
}

/* ================== SmithWatermanFitness ================== */
double SmithWatermanFitness::calculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!individual || !representation) {
        LOG_ERROR("Null individual or representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    try {
        std::string dnaStr = representation->toDNA(individual, instance);
        if (dnaStr.empty()) {
            LOG_ERROR("Empty DNA sequence generated from individual");
            return -std::numeric_limits<double>::infinity();
        }

        const std::string& original = instance.getOriginalDNA();
        return static_cast<double>(smithWaterman(dnaStr, original));
    } catch (const std::exception& e) {
        LOG_ERROR("Error calculating Smith-Waterman fitness: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

int SmithWatermanFitness::smithWaterman(const std::string& seq1, const std::string& seq2) const
{
    if (seq1.empty() || seq2.empty()) return 0;

    static const int MATCH_SCORE = 2;
    static const int MISMATCH_SCORE = -1;
    static const int GAP_PENALTY = -2;

    std::vector<std::vector<int>> dp(seq1.length() + 1, 
                                   std::vector<int>(seq2.length() + 1, 0));

    int maxScore = 0;

    for (size_t i = 1; i <= seq1.length(); i++) {
        for (size_t j = 1; j <= seq2.length(); j++) {
            const int match = dp[i-1][j-1] + 
                (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            const int del = dp[i-1][j] + GAP_PENALTY;
            const int ins = dp[i][j-1] + GAP_PENALTY;

            dp[i][j] = std::max({0, match, del, ins});
            maxScore = std::max(maxScore, dp[i][j]);
        }
    }

    return maxScore;
}

/* ================== OptimizedGraphBasedFitness ================== */
double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!individual || !representation) {
        LOG_ERROR("Null individual or representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    try {
        const auto& genes = individual->getGenes();
        if (genes.empty()) return 0.0;

        std::lock_guard<std::mutex> lock(graphCacheMutex);
        ensureBufferSizes(genes.size());
        buildGraph(genes, instance);

        const double connectivityScore = calculateConnectivityScore();
        const double coverageScore = calculateSpectrumCoverageScore(genes, instance);
        const double lengthPenalty = calculateLengthPenalty(genes, instance);

        return connectivityScore * coverageScore * (1.0 - lengthPenalty);
    } catch (const std::exception& e) {
        LOG_ERROR("Error calculating graph-based fitness: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

void OptimizedGraphBasedFitness::ensureBufferSizes(size_t size) const
{
    adjacencyMatrix.resize(size, std::vector<bool>(size, false));
    distanceMatrix.resize(size, std::vector<int>(size, 0));
}

void OptimizedGraphBasedFitness::buildGraph(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    const size_t n = genes.size();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            adjacencyMatrix[i][j] = areNodesConnected(genes[i], genes[j], instance);
            distanceMatrix[i][j] = adjacencyMatrix[i][j] ? 1 : std::numeric_limits<int>::max();
        }
    }

    // Floyd-Warshall algorithm for all-pairs shortest paths
    for (size_t k = 0; k < n; ++k) {
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                if (distanceMatrix[i][k] != std::numeric_limits<int>::max() && 
                    distanceMatrix[k][j] != std::numeric_limits<int>::max()) {
                    distanceMatrix[i][j] = std::min(distanceMatrix[i][j], 
                                                  distanceMatrix[i][k] + distanceMatrix[k][j]);
                }
            }
        }
    }
}

bool OptimizedGraphBasedFitness::areNodesConnected(
    int node1,
    int node2,
    const DNAInstance& instance) const
{
    if (node1 < 0 || node2 < 0) return false;
    
    const auto& spectrum = instance.getSpectrum();
    if (node1 >= static_cast<int>(spectrum.size()) || 
        node2 >= static_cast<int>(spectrum.size())) {
        return false;
    }

    const std::string& kmer1 = spectrum[node1];
    const std::string& kmer2 = spectrum[node2];
    const int k = instance.getK();

    if (static_cast<int>(kmer1.length()) < k - 1 || 
        static_cast<int>(kmer2.length()) < k - 1) {
        return false;
    }

    std::string suffix = kmer1.substr(kmer1.length() - (k - 1));
    std::string prefix = kmer2.substr(0, k - 1);

    return suffix == prefix;
}

double OptimizedGraphBasedFitness::calculateConnectivityScore() const
{
    const size_t n = adjacencyMatrix.size();
    if (n == 0) return 0.0;

    double totalConnectivity = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j && distanceMatrix[i][j] != std::numeric_limits<int>::max()) {
                totalConnectivity += 1.0 / distanceMatrix[i][j];
            }
        }
    }

    return totalConnectivity / (n * (n - 1));
}

double OptimizedGraphBasedFitness::calculateSpectrumCoverageScore(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    const auto& spectrum = instance.getSpectrum();
    std::vector<bool> covered(spectrum.size(), false);
    double coverage = 0.0;

    for (int gene : genes) {
        if (gene >= 0 && gene < static_cast<int>(spectrum.size()) && !covered[gene]) {
            covered[gene] = true;
            coverage += 1.0;
        }
    }

    return coverage / spectrum.size();
}

double OptimizedGraphBasedFitness::calculateLengthPenalty(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    const int targetLength = instance.getN();
    if (targetLength <= 0) return 0.0;

    const auto& spectrum = instance.getSpectrum();
    int actualLength = 0;
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] >= 0 && genes[i] < static_cast<int>(spectrum.size())) {
            actualLength += spectrum[genes[i]].length();
            if (i > 0) {
                actualLength -= instance.getK() - 1;
            }
        }
    }

    double lengthDiff = std::abs(actualLength - targetLength);
    return std::exp(-lengthDiff / targetLength);
}

void OptimizedGraphBasedFitness::cleanupCache() const
{
    // Clear matrices when they get too large
    if (adjacencyMatrix.size() > 1000) {
        adjacencyMatrix.clear();
        distanceMatrix.clear();
    }
}
