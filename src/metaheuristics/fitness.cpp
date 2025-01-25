//
// Created by konrad_guest on 07/01/2025.
// SMART

#include <unordered_set>
#include <iostream>
#include "../include/metaheuristics/fitness.h"
#include "../include/utils/logging.h"
#include <algorithm>
#include <numeric>
#include <sstream>
#include "../../include/metaheuristics/fitness_impl.h"
#include <cmath>

/* ================== SimpleFitness ================== */
double SimpleFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    // Convert solution to DNA sequence
    std::vector<char> dna = representation->toDNA(solution, instance);
    if (dna.empty()) {
        LOG_ERROR("Empty DNA sequence generated from solution");
        return -std::numeric_limits<double>::infinity();
    }

    // Calculate fitness components
    double connectivityScore = calculateConnectivityScore(solution, instance);
    double spectrumCoverageScore = calculateSpectrumCoverageScore(dna, instance);
    double lengthPenalty = calculateLengthPenalty(dna.size(), instance.getN());

    // Combine scores with weights
    constexpr double CONNECTIVITY_WEIGHT = 0.4;
    constexpr double COVERAGE_WEIGHT = 0.4;
    constexpr double LENGTH_WEIGHT = 0.2;

    double totalFitness = (CONNECTIVITY_WEIGHT * connectivityScore) +
                         (COVERAGE_WEIGHT * spectrumCoverageScore) +
                         (LENGTH_WEIGHT * lengthPenalty);

    return totalFitness;
}

double SimpleFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;

    double score = 0.0;
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();

    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= (int)spectrum.size() ||
            genes[i + 1] < 0 || genes[i + 1] >= (int)spectrum.size()) {
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i + 1]];

        if (current.length() >= k - 1 && next.length() >= k - 1) {
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
    int k = instance.getK();
    
    if (dna.size() < k) return 0.0;

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
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    // Convert solution to DNA sequence
    std::vector<char> dna = representation->toDNA(solution, instance);
    if (dna.empty()) {
        LOG_ERROR("Empty DNA sequence generated from solution");
        return -std::numeric_limits<double>::infinity();
    }

    // Calculate fitness components
    double connectivityScore = calculateConnectivityScore(solution, instance);
    double spectrumCoverageScore = calculateSpectrumCoverageScore(dna, instance);
    double lengthPenalty = calculateLengthPenalty(dna.size(), instance.getN());
    double sequenceSimilarity = calculateSequenceSimilarity(dna, instance);

    // Combine scores with weights
    constexpr double CONNECTIVITY_WEIGHT = 0.35;
    constexpr double COVERAGE_WEIGHT = 0.35;
    constexpr double LENGTH_WEIGHT = 0.15;
    constexpr double SIMILARITY_WEIGHT = 0.15;

    double totalFitness = (CONNECTIVITY_WEIGHT * connectivityScore) +
                         (COVERAGE_WEIGHT * spectrumCoverageScore) +
                         (LENGTH_WEIGHT * lengthPenalty) +
                         (SIMILARITY_WEIGHT * sequenceSimilarity);

    return totalFitness;
}

double BetterFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;

    double score = 0.0;
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();

    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= (int)spectrum.size() ||
            genes[i + 1] < 0 || genes[i + 1] >= (int)spectrum.size()) {
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i + 1]];

        if (current.length() >= k - 1 && next.length() >= k - 1) {
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
    int k = instance.getK();
    
    if (dna.size() < k) return 0.0;

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

    // Convert dna vector to string
    std::string dnaStr(dna.begin(), dna.end());

    // Calculate Levenshtein distance
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
            if (dnaStr[i-1] == originalDNA[j-1]) {
                dp[i][j] = dp[i-1][j-1];
            } else {
                dp[i][j] = 1 + std::min({dp[i-1][j],      // deletion
                                       dp[i][j-1],         // insertion
                                       dp[i-1][j-1]});     // substitution
            }
        }
    }

    int distance = dp[dnaStr.length()][originalDNA.length()];
    double maxDistance = std::max(dnaStr.length(), originalDNA.length());
    
    // Convert distance to similarity score (0 to 1)
    return 1.0 - (distance / maxDistance);
}

/* ================== SmithWatermanFitness ================== */
double SmithWatermanFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    // Convert solution to DNA sequence
    std::vector<char> dna = representation->toDNA(solution, instance);
    if (dna.empty()) {
        LOG_ERROR("Empty DNA sequence generated from solution");
        return -std::numeric_limits<double>::infinity();
    }

    // Convert DNA vector to string
    std::string dnaStr(dna.begin(), dna.end());
    
    // Calculate k-mer matches
    double kmerScore = 0.0;
    const auto& spectrum = instance.getSpectrum();
    for (const auto& kmer : spectrum) {
        if (dnaStr.find(kmer) != std::string::npos) {
            kmerScore += 1.0;
        }
    }
    kmerScore /= spectrum.size();

    // Calculate Smith-Waterman similarity to original DNA
    const std::string& originalDNA = instance.getOriginalDNA();
    int similarity = smithWaterman(dnaStr, originalDNA);
    double similarityScore = static_cast<double>(similarity) / 
                           std::max(dnaStr.length(), originalDNA.length());

    // Combine scores with weights
    constexpr double KMER_WEIGHT = 0.7;
    constexpr double SIMILARITY_WEIGHT = 0.3;
    
    return (KMER_WEIGHT * kmerScore) + (SIMILARITY_WEIGHT * similarityScore);
}

int SmithWatermanFitness::smithWaterman(const std::string& seq1, const std::string& seq2) const {
    std::vector<std::vector<int>> matrix(seq1.length() + 1, std::vector<int>(seq2.length() + 1, 0));
    int maxScore = 0;

    for (size_t i = 1; i <= seq1.length(); ++i) {
        for (size_t j = 1; j <= seq2.length(); ++j) {
            int match = matrix[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            int del = matrix[i-1][j] + GAP_PENALTY;
            int ins = matrix[i][j-1] + GAP_PENALTY;
            
            matrix[i][j] = std::max({0, match, del, ins});
            maxScore = std::max(maxScore, matrix[i][j]);
        }
    }

    return maxScore;
}

/* ================== OptimizedGraphBasedFitness ================== */
double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    // Convert solution to DNA sequence
    std::vector<char> dna = representation->toDNA(solution, instance);
    if (dna.empty()) {
        LOG_ERROR("Empty DNA sequence generated from solution");
        return -std::numeric_limits<double>::infinity();
    }

    // Calculate fitness components
    double connectivityScore = calculateConnectivityScore(solution, instance);
    double spectrumCoverageScore = calculateSpectrumCoverageScore(dna, instance);
    double lengthPenalty = calculateLengthPenalty(dna.size(), instance.getN());
    double pathOptimalityScore = calculatePathOptimalityScore(solution, instance);

    // Combine scores with weights
    constexpr double CONNECTIVITY_WEIGHT = 0.35;
    constexpr double COVERAGE_WEIGHT = 0.35;
    constexpr double LENGTH_WEIGHT = 0.15;
    constexpr double OPTIMALITY_WEIGHT = 0.15;

    double totalFitness = (CONNECTIVITY_WEIGHT * connectivityScore) +
                         (COVERAGE_WEIGHT * spectrumCoverageScore) +
                         (LENGTH_WEIGHT * lengthPenalty) +
                         (OPTIMALITY_WEIGHT * pathOptimalityScore);

    return totalFitness;
}

double OptimizedGraphBasedFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;

    double score = 0.0;
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();

    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= (int)spectrum.size() ||
            genes[i + 1] < 0 || genes[i + 1] >= (int)spectrum.size()) {
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i + 1]];

        if (current.length() >= k - 1 && next.length() >= k - 1) {
            std::string suffix = current.substr(current.length() - (k - 1));
            std::string prefix = next.substr(0, k - 1);
            if (suffix == prefix) {
                score += 1.0;
            }
        }
    }

    return score;
}

double OptimizedGraphBasedFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const
{
    if (dna.empty()) return 0.0;

    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    if (dna.size() < k) return 0.0;

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

double OptimizedGraphBasedFitness::calculateLengthPenalty(int actualLength, int targetLength) const
{
    if (targetLength <= 0) return 0.0;
    
    double lengthDiff = std::abs(actualLength - targetLength);
    return std::exp(-lengthDiff / targetLength);
}

double OptimizedGraphBasedFitness::calculatePathOptimalityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;

    // Build adjacency matrix for the graph
    auto adjacencyMatrix = buildAdjacencyMatrix(instance);
    
    // Calculate path score
    double pathScore = 0.0;
    double maxPossibleScore = 0.0;
    
    // Calculate actual path score
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= (int)adjacencyMatrix.size() ||
            genes[i + 1] < 0 || genes[i + 1] >= (int)adjacencyMatrix.size()) {
            continue;
        }
        
        const auto& edge = adjacencyMatrix[genes[i]][genes[i + 1]];
        if (edge.exists) {
            pathScore += edge.weight;
        }
    }
    
    // Calculate maximum possible score for path of this length
    for (const auto& row : adjacencyMatrix) {
        double maxEdgeWeight = 0.0;
        for (const auto& edge : row) {
            if (edge.exists) {
                maxEdgeWeight = std::max(maxEdgeWeight, edge.weight);
            }
        }
        maxPossibleScore += maxEdgeWeight;
    }
    
    return maxPossibleScore > 0.0 ? pathScore / maxPossibleScore : 0.0;
}

std::vector<std::vector<PreprocessedEdge>> OptimizedGraphBasedFitness::buildAdjacencyMatrix(
    const DNAInstance& instance) const
{
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    std::vector<std::vector<PreprocessedEdge>> matrix(
        spectrum.size(),
        std::vector<PreprocessedEdge>(spectrum.size())
    );
    
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
                if (weight > 0) {
                    matrix[i][j] = PreprocessedEdge(j, weight, true);
                }
            }
        }
    }
    
    return matrix;
}

int OptimizedGraphBasedFitness::calculateEdgeWeight(
    const std::string& from,
    const std::string& to,
    int k) const
{
    if (from.length() < k - 1 || to.length() < k - 1) return 0;
    
    // Compare suffix of 'from' with prefix of 'to'
    std::string suffix = from.substr(from.length() - (k - 1));
    std::string prefix = to.substr(0, k - 1);
    
    if (suffix == prefix) {
        return k - 1;  // Weight is the length of the overlap
    }
    
    return 0;
}

double Fitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to calculateFitness");
        return -std::numeric_limits<double>::infinity();
    }

    // Check cache first if available
    if (m_cache) {
        return m_cache->getOrCalculateFitness(solution, instance, shared_from_this(), representation);
    }

    // Convert to DNA and calculate fitness
    std::vector<char> dna = representation->toDNA(solution, instance);
    if (dna.empty()) {
        LOG_ERROR("Empty DNA sequence generated from solution");
        return -std::numeric_limits<double>::infinity();
    }

    return calculateDNAFitness(dna, instance);
}

double Fitness::calculateDNAFitness(const std::vector<char>& dna, const DNAInstance& instance) const {
    if (dna.empty()) {
        LOG_ERROR("Empty DNA sequence provided to calculateDNAFitness");
        return -std::numeric_limits<double>::infinity();
    }

    // Calculate Levenshtein distance between reconstructed DNA and original DNA
    const std::string& originalDNA = instance.getOriginalDNA();
    const std::string reconstructedDNA(dna.begin(), dna.end());
    
    int distance = levenshteinDistance(originalDNA, reconstructedDNA);
    
    // Convert distance to fitness (lower distance = higher fitness)
    // We use negative distance so that higher values are better
    return -static_cast<double>(distance);
}

int Fitness::levenshteinDistance(const std::string& s1, const std::string& s2) const {
    if (s1.empty()) return s2.length();
    if (s2.empty()) return s1.length();

    const std::size_t len1 = s1.size(), len2 = s2.size();
    std::vector<std::vector<int>> d(len1 + 1, std::vector<int>(len2 + 1));

    for (int i = 0; i <= len1; ++i)
        d[i][0] = i;
    for (int j = 0; j <= len2; ++j)
        d[0][j] = j;

    for (int i = 1; i <= len1; ++i)
        for (int j = 1; j <= len2; ++j)
            d[i][j] = std::min({ d[i - 1][j] + 1,
                                d[i][j - 1] + 1,
                                d[i - 1][j - 1] + (s1[i - 1] != s2[j - 1]) });

    return d[len1][len2];
}

double OptimizedGraphBasedFitness::evaluate(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance)
{
    if (!individual) {
        LOG_ERROR("Cannot evaluate null individual");
        return -std::numeric_limits<double>::infinity();
    }
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) {
        LOG_ERROR("Cannot evaluate empty individual");
        return -std::numeric_limits<double>::infinity();
    }
    
    try {
        // Calculate hash for caching
        std::string key = individual->toString();
        
        // Check cache first
        {
            std::lock_guard<std::mutex> lock(graphCacheMutex);
            auto it = cache.find(key);
            if (it != cache.end()) {
                return it->second;
            }
        }
        
        // Ensure buffers are properly sized
        ensureBufferSizes(instance.getSpectrum().size());
        
        // Build graph using pre-allocated buffers
        buildGraph(genes, instance);
        
        // Calculate fitness components
        double connectivityScore = calculateConnectivityScore();
        double spectrumCoverageScore = calculateSpectrumCoverageScore(genes, instance);
        double lengthPenalty = calculateLengthPenalty(genes, instance);
        
        // Combine scores with weights
        double fitness = (0.5 * connectivityScore + 
                        0.3 * spectrumCoverageScore - 
                        0.2 * lengthPenalty);
        
        // Cache result with thread safety
        {
            std::lock_guard<std::mutex> lock(graphCacheMutex);
            if (cache.size() >= maxCacheSize) {
                cleanupCache();
            }
            cache[key] = fitness;
        }
        
        return fitness;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Fitness evaluation failed: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

void OptimizedGraphBasedFitness::ensureBufferSizes(size_t size) {
    if (nodeUsageBuffer.size() < size) {
        nodeUsageBuffer.resize(size, 0);
    }
    if (coveredBuffer.size() < size) {
        coveredBuffer.resize(size, false);
    }
    if (adjacencyMatrixBuffer.size() < size) {
        adjacencyMatrixBuffer.resize(size, std::vector<PreprocessedEdge>(size));
    }
    
    // Clear buffers
    std::fill(nodeUsageBuffer.begin(), nodeUsageBuffer.end(), 0);
    std::fill(coveredBuffer.begin(), coveredBuffer.end(), false);
    for (auto& row : adjacencyMatrixBuffer) {
        std::fill(row.begin(), row.end(), PreprocessedEdge{});
    }
}

void OptimizedGraphBasedFitness::buildGraph(
    const std::vector<int>& genes,
    const DNAInstance& instance)
{
    const auto& spectrum = instance.getSpectrum();
    size_t n = genes.size();
    
    for (size_t i = 0; i < n; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) {
            continue;
        }
        for (size_t j = i + 1; j < n; ++j) {
            if (genes[j] < 0 || genes[j] >= static_cast<int>(spectrum.size())) {
                continue;
            }
            if (areNodesConnected(genes[i], genes[j], instance)) {
                adjacencyMatrixBuffer[i][j].exists = true;
                adjacencyMatrixBuffer[i][j].weight = 1;
                adjacencyMatrixBuffer[j][i].exists = true;
                adjacencyMatrixBuffer[j][i].weight = 1;
            }
        }
    }
}

bool OptimizedGraphBasedFitness::areNodesConnected(
    int node1,
    int node2,
    const DNAInstance& instance) const
{
    const auto& spectrum = instance.getSpectrum();
    if (node1 < 0 || node1 >= static_cast<int>(spectrum.size()) ||
        node2 < 0 || node2 >= static_cast<int>(spectrum.size())) {
        return false;
    }

    const std::string& kmer1 = spectrum[node1];
    const std::string& kmer2 = spectrum[node2];
    
    if (kmer1.length() < instance.getK() - 1 || kmer2.length() < instance.getK() - 1) {
        return false;
    }

    std::string suffix = kmer1.substr(kmer1.length() - (instance.getK() - 1));
    std::string prefix = kmer2.substr(0, instance.getK() - 1);
    
    return suffix == prefix;
}

double OptimizedGraphBasedFitness::calculateConnectivityScore() const {
    if (adjacencyMatrixBuffer.empty()) {
        return 0.0;
    }
    
    size_t totalEdges = 0;
    size_t maxPossibleEdges = (adjacencyMatrixBuffer.size() * 
                              (adjacencyMatrixBuffer.size() - 1)) / 2;
    
    for (size_t i = 0; i < adjacencyMatrixBuffer.size(); ++i) {
        for (size_t j = i + 1; j < adjacencyMatrixBuffer.size(); ++j) {
            if (adjacencyMatrixBuffer[i][j].exists) {
                ++totalEdges;
            }
        }
    }
    
    return maxPossibleEdges > 0 ? static_cast<double>(totalEdges) / maxPossibleEdges : 0.0;
}

double OptimizedGraphBasedFitness::calculateSpectrumCoverageScore(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    if (genes.empty()) {
        return 0.0;
    }
    
    const auto& spectrum = instance.getSpectrum();
    size_t coveredCount = 0;
    
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) {
            continue;
        }
        if (!coveredBuffer[genes[i]]) {
            coveredBuffer[genes[i]] = true;
            ++coveredCount;
        }
    }
    
    return spectrum.size() > 0 ? static_cast<double>(coveredCount) / spectrum.size() : 0.0;
}

double OptimizedGraphBasedFitness::calculateLengthPenalty(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    double expectedLength = static_cast<double>(instance.getN());
    double actualLength = static_cast<double>(genes.size());
    double lengthDiff = std::abs(actualLength - expectedLength);
    
    // Normalize penalty to [0,1] range
    return std::exp(-lengthDiff / expectedLength);
}

void OptimizedGraphBasedFitness::cleanupCache() {
    if (cache.size() <= maxCacheSize / 2) return;
    
    std::vector<std::pair<std::string, double>> entries(
        cache.begin(), cache.end());
    
    // Keep only the most recent half
    size_t keepCount = maxCacheSize / 2;
    cache.clear();
    for (size_t i = entries.size() - keepCount; i < entries.size(); ++i) {
        cache.insert(entries[i]);
    }
}
