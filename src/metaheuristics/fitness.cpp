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

/* ================== SimpleFitness ================== */
double SimpleFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!solution || solution->empty()) {
        LOG_WARNING("Null or empty solution provided to fitness calculation");
        return 0.0;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to fitness calculation");
        return 0.0;
    }

    // Convert solution to DNA string
    std::string dna = representation->toString(solution, instance);
    
    // Count matching k-mers
    int matches = 0;
    const auto& spectrum = instance.getSpectrum();
    for (const auto& kmer : spectrum) {
        if (dna.find(kmer) != std::string::npos) {
            matches++;
        }
    }

    return static_cast<double>(matches) / spectrum.size();
}

/* ================== BetterFitness ================== */
double BetterFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!solution || solution->empty()) {
        LOG_WARNING("Null or empty solution provided to fitness calculation");
        return 0.0;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to fitness calculation");
        return 0.0;
    }

    // Get DNA string representation
    std::string dna = representation->toString(solution, instance);
    
    // Calculate k-mer matches
    double kmerScore = 0.0;
    const auto& spectrum = instance.getSpectrum();
    for (const auto& kmer : spectrum) {
        if (dna.find(kmer) != std::string::npos) {
            kmerScore += 1.0;
        }
    }
    kmerScore /= spectrum.size();

    // Calculate Levenshtein distance to original DNA
    const std::string& originalDNA = instance.getOriginalDNA();
    int distance = levenshteinDistance(dna, originalDNA);
    double distanceScore = 1.0 - (static_cast<double>(distance) / std::max(dna.length(), originalDNA.length()));

    // Combine scores with weights
    constexpr double KMER_WEIGHT = 0.7;
    constexpr double DISTANCE_WEIGHT = 0.3;
    return (KMER_WEIGHT * kmerScore) + (DISTANCE_WEIGHT * distanceScore);
}

/* ================== SmithWatermanFitness ================== */
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

double SmithWatermanFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!solution || solution->empty()) {
        LOG_WARNING("Null or empty solution provided to fitness calculation");
        return 0.0;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to fitness calculation");
        return 0.0;
    }

    // Get DNA string representation
    std::string dna = representation->toString(solution, instance);
    
    // Calculate Smith-Waterman score
    int swScore = smithWaterman(dna, instance.getOriginalDNA());
    
    // Normalize score
    double maxPossibleScore = std::min(dna.length(), instance.getOriginalDNA().length()) * MATCH_SCORE;
    return static_cast<double>(swScore) / maxPossibleScore;
}

/* ================== OptimizedGraphBasedFitness ================== */

// Pomocnicza inicjalizacja bufora
void OptimizedGraphBasedFitness::initBuffers(size_t size) const {
    if (nodeUsageBuffer.size() < size) {
        nodeUsageBuffer.resize(size);
    }
}

// Generowanie klucza do kesza
std::string OptimizedGraphBasedFitness::createCacheKey(const std::vector<std::string>& spectrum, int k) const {
    std::string key = std::to_string(k) + ":";
    for (const auto& s : spectrum) {
        key += s + ",";
    }
    return key;
}

// Właściwe budowanie grafu (lub użycie z kesza)
std::vector<std::vector<OptimizedGraphBasedFitness::Edge>> 
OptimizedGraphBasedFitness::buildSpectrumGraph(
    const std::vector<std::string>& spectrum, 
    int k) const 
{
    std::string cacheKey = createCacheKey(spectrum, k);
    auto it = graphCache.find(cacheKey);
    if (it != graphCache.end()) {
        return it->second;
    }
    
    int n = (int)spectrum.size();
    std::vector<std::vector<Edge>> graph(n);
    
    // Równoległe wypełnianie list sąsiedztwa
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        graph[i].reserve(n/4); // np. heurystycznie
        for (int j = 0; j < n; j++) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
                if (weight > 0) {
                    // W sekcji równoległej dostęp do graph[i] jest bezpieczny,
                    // bo i jest unikalne w każdej iteracji i.
                    graph[i].emplace_back(j, weight);
                }
            }
        }
    }
    
    // Zapamiętujemy w keszu
    graphCache[cacheKey] = graph;
    return graph;
}

// Wyznaczenie wagi krawędzi między dwoma k-merami
int OptimizedGraphBasedFitness::calculateEdgeWeight(
    const std::string& from, 
    const std::string& to, 
    int k) const 
{
    if (from.length() < k || to.length() < k) {
        return 0;
    }

    // Sprawdzamy nakładanie się końca pierwszego k-meru z początkiem drugiego
    int maxOverlap = 0;
    for (int overlap = k-1; overlap > 0; overlap--) {
        if (from.substr(from.length()-overlap) == to.substr(0, overlap)) {
            maxOverlap = overlap;
            break;
        }
    }

    // Waga krawędzi zależy od długości nakładania
    if (maxOverlap == k-1) return 3;  // Idealne nakładanie
    if (maxOverlap >= k/2) return 2;  // Dobre nakładanie
    if (maxOverlap > 0) return 1;     // Słabe nakładanie
    return 0;                         // Brak nakładania
}

// Analiza ścieżki w grafie: w tej nowej wersji adjacencyMatrix budujemy lokalnie.
OptimizedGraphBasedFitness::PathAnalysis 
OptimizedGraphBasedFitness::analyzePath(
    const std::vector<int>& path,
    const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix
) const {
    PathAnalysis analysis(0, 0, 0, 0);
    std::fill(nodeUsageBuffer.begin(), nodeUsageBuffer.end(), 0);

    for (size_t i = 0; i < path.size(); ++i) {
        int node = path[i];
        nodeUsageBuffer[node]++;

        if (i > 0) {
            int prevNode = path[i-1];
            const auto& edge = adjacencyMatrix[prevNode][node];
            
            if (edge.exists) {
                if (edge.weight == 1) {
                    analysis.edgesWeight1++;
                } else {
                    analysis.edgesWeight2or3++;
                }
            }
        }
    }

    for (int usage : nodeUsageBuffer) {
        if (usage > 0) {
            analysis.uniqueNodesUsed++;
            if (usage > 1) {
                analysis.repeatNodeUsages += usage - 1;
            }
        }
    }

    return analysis;
}

// Prosta funkcja – w permutacji 1:1 to nasza ścieżka.
const std::vector<int>& 
OptimizedGraphBasedFitness::permutationToPath(std::shared_ptr<std::vector<int>> individual) const {
    return *individual;
}

double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!solution || solution->empty()) {
        LOG_WARNING("Null or empty solution provided to fitness calculation");
        return 0.0;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to fitness calculation");
        return 0.0;
    }

    // Initialize buffers if needed
    initBuffers(instance.getSpectrum().size());

    // Build or get cached graph
    std::string cacheKey = createCacheKey(instance.getSpectrum(), instance.getK());
    auto& graph = graphCache[cacheKey];
    if (graph.empty()) {
        graph = buildSpectrumGraph(instance.getSpectrum(), instance.getK());
    }

    // Convert graph to adjacency matrix for faster access
    auto adjacencyMatrix = buildAdjacencyMatrix(graph);

    // Analyze path properties
    const auto& path = permutationToPath(solution);
    auto analysis = analyzePath(path, adjacencyMatrix);

    // Calculate fitness components
    double connectivityScore = (analysis.edgesWeight1 + 0.5 * analysis.edgesWeight2or3) / 
                             static_cast<double>(path.size() - 1);
    double uniquenessScore = static_cast<double>(analysis.uniqueNodesUsed) / instance.getSpectrum().size();
    double repetitionPenalty = 1.0 / (1.0 + analysis.repeatNodeUsages);

    // Combine scores with weights
    constexpr double CONNECTIVITY_WEIGHT = 0.5;
    constexpr double UNIQUENESS_WEIGHT = 0.3;
    constexpr double REPETITION_WEIGHT = 0.2;

    return (CONNECTIVITY_WEIGHT * connectivityScore) +
           (UNIQUENESS_WEIGHT * uniquenessScore) +
           (REPETITION_WEIGHT * repetitionPenalty);
}

std::vector<std::vector<OptimizedGraphBasedFitness::PreprocessedEdge>> OptimizedGraphBasedFitness::buildAdjacencyMatrix(
    const std::vector<std::vector<Edge>>& graph
) const {
    std::vector<std::vector<PreprocessedEdge>> matrix(
        graph.size(),
        std::vector<PreprocessedEdge>(graph.size())
    );

    for (size_t i = 0; i < graph.size(); ++i) {
        for (const auto& edge : graph[i]) {
            matrix[i][edge.to] = PreprocessedEdge(edge.to, edge.weight, true);
        }
    }

    return matrix;
}

double Fitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!individual || individual->empty()) {
        LOG_WARNING("Null or empty individual provided to fitness calculation");
        return 0.0;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to fitness calculation");
        return 0.0;
    }

    // Check cache first if available
    if (m_cache) {
        return m_cache->getOrCalculateFitness(individual, instance, shared_from_this(), representation);
    }

    // Convert to DNA and calculate fitness
    std::vector<char> dna = representation->toDNA(individual, instance);
    return calculateDNAFitness(dna, instance);
}

double Fitness::calculateDNAFitness(const std::vector<char>& dna, const DNAInstance& instance) const {
    // Calculate Levenshtein distance between reconstructed DNA and original DNA
    const std::string& originalDNA = instance.getDNA();
    const std::string reconstructedDNA(dna.begin(), dna.end());
    
    int distance = levenshteinDistance(originalDNA, reconstructedDNA);
    
    // Convert distance to fitness (lower distance = higher fitness)
    // We use negative distance so that higher values are better
    return -static_cast<double>(distance);
}

int Fitness::levenshteinDistance(const std::string& s1, const std::string& s2) const {
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
