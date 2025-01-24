//
// Created by konrad_guest on 07/01/2025.
// SMART

#include <unordered_set>
#include <iostream>
#include "metaheuristics/fitness.h"
#include "utils/logging.h"
#include <algorithm>
#include <numeric>
#include <sstream>

/* ================== SimpleFitness ================== */
double SimpleFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!individual || individual->empty()) {
        return 0.0;
    }

    // Convert individual to DNA string
    std::string dna;
    for (int gene : *individual) {
        if (gene >= 0 && gene < static_cast<int>(instance.getSpectrum().size())) {
            dna += instance.getSpectrum()[gene];
        } else {
            LOG_ERROR("Invalid gene value: " + std::to_string(gene));
            return 0.0;
        }
    }

    // Count matching k-mers
    int k = instance.getK();
    int matches = 0;
    const auto& spectrum = instance.getSpectrum();
    
    for (size_t i = 0; i <= dna.length() - k; ++i) {
        std::string kmer = dna.substr(i, k);
        if (std::find(spectrum.begin(), spectrum.end(), kmer) != spectrum.end()) {
            matches++;
        }
    }

    return static_cast<double>(matches);
}

/* ================== BetterFitness ================== */
double BetterFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!individual || individual->empty()) {
        return 0.0;
    }

    // Get base fitness from k-mer matches
    double baseFitness = SimpleFitness().calculateFitness(individual, instance, representation);
    
    // Add bonus for sequence length close to expected
    int expectedLength = instance.getSize();
    int actualLength = individual->size() * instance.getK();
    double lengthPenalty = 1.0 - std::abs(expectedLength - actualLength) / static_cast<double>(expectedLength);
    
    return baseFitness * (0.8 + 0.2 * lengthPenalty);
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
    const std::shared_ptr<std::vector<int>>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!individual || individual->empty()) {
        return 0.0;
    }

    // Convert individual to DNA string
    std::string dna;
    for (int gene : *individual) {
        if (gene >= 0 && gene < static_cast<int>(instance.getSpectrum().size())) {
            dna += instance.getSpectrum()[gene];
        } else {
            LOG_ERROR("Invalid gene value: " + std::to_string(gene));
            return 0.0;
        }
    }

    // Calculate Smith-Waterman score against each k-mer
    double totalScore = 0.0;
    const auto& spectrum = instance.getSpectrum();
    
    for (const auto& kmer : spectrum) {
        totalScore += smithWaterman(dna, kmer);
    }

    return totalScore;
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
    // Podobna logika co wcześniej
    const char* fromEnd = from.data() + from.size() - (k - 1);
    const char* toStart = to.data();
    
    bool match = true;
    for (int i = 0; i < k - 1; i++) {
        if (fromEnd[i] != toStart[i]) {
            match = false;
            break;
        }
    }
    if (match) return 1;
    
    // sprawdzamy jeszcze "kolejne" przesunięcia ...
    fromEnd++;
    match = true;
    for (int i = 0; i < k - 2; i++) {
        if (fromEnd[i] != toStart[i]) {
            match = false;
            break;
        }
    }
    if (match) return 2;
    
    fromEnd++;
    // Jeszcze minimalna zbieżność – 1 znak
    if (fromEnd[0] == toStart[0]) return 3;
    
    return 0;
}

// Analiza ścieżki w grafie: w tej nowej wersji adjacencyMatrix budujemy lokalnie.
OptimizedGraphBasedFitness::PathAnalysis 
OptimizedGraphBasedFitness::analyzePath(const std::vector<int>& path,
                                        const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix) const 
{
    // Wyzerowanie bufora
    std::fill(nodeUsageBuffer.begin(), nodeUsageBuffer.end(), 0);

    PathAnalysis analysis(nodeUsageBuffer);

    for (size_t i = 0; i < path.size() - 1; i++) {
        int from = path[i];
        int to = path[i + 1];

        // **Debug Log**: Przed dostępem do adjacencyMatrix
        std::cout << "[DEBUG analyzePath] Edge from " << from << " to " << to << "\n";

        if (from < 0 || from >= (int)adjacencyMatrix.size() ||
            to   < 0 || to   >= (int)adjacencyMatrix.size())
        {
            std::cerr << "[ERROR analyzePath] Edge indices out-of-bounds: from=" << from << ", to=" << to << "\n";
            continue; // lub inna logika obsługi błędu
        }
        
        if (from >= 0 && from < (int)adjacencyMatrix.size() &&
            to   >= 0 && to   < (int)adjacencyMatrix.size())
        {
            // zliczamy użycie węzła
            analysis.nodeUsageCount[from]++;
            if (analysis.nodeUsageCount[from] == 1) {
                analysis.uniqueNodesUsed++;
            } else {
                analysis.repeatNodeUsages++;
            }

            // sprawdzamy krawędź
            const auto& edge = adjacencyMatrix[from][to];
            if (edge.exists) {
                if (edge.weight == 1) {
                    analysis.edgesWeight1++;
                } else {
                    analysis.edgesWeight2or3++;
                }
            }
        }
    }
    // Ostatni węzeł
    int last = path.back();
    if (!path.empty()) {
        if (last >= 0 && last < (int)adjacencyMatrix.size()) {
            analysis.nodeUsageCount[last]++;
            if (analysis.nodeUsageCount[last] == 1) {
                analysis.uniqueNodesUsed++;
            } else {
                analysis.repeatNodeUsages++;
            }
        }
    }else
    {
        std::cerr << "[ERROR analyzePath] Last node out-of-bounds: " << last << "\n";
    }
    // **Debug Log**: Wyniki analizy ścieżki
    std::cout << "[DEBUG analyzePath] Analysis results - edgesWeight1: " << analysis.edgesWeight1
              << ", edgesWeight2or3: " << analysis.edgesWeight2or3
              << ", uniqueNodesUsed: " << analysis.uniqueNodesUsed
              << ", repeatNodeUsages: " << analysis.repeatNodeUsages << "\n";

    return analysis;
}

// Prosta funkcja – w permutacji 1:1 to nasza ścieżka.
const std::vector<int>& 
OptimizedGraphBasedFitness::permutationToPath(std::shared_ptr<std::vector<int>> individual) const {
    return *individual;
}

double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!individual || individual->empty()) {
        return 0.0;
    }

    // Build spectrum graph
    auto graph = buildSpectrumGraph(instance.getSpectrum(), instance.getK());
    auto adjacencyMatrix = buildAdjacencyMatrix(graph);
    
    // Analyze path
    const auto& path = permutationToPath(individual);
    auto analysis = analyzePath(path, adjacencyMatrix);
    
    // Calculate fitness based on path analysis
    double pathScore = analysis.edgesWeight1 * 1.0 + analysis.edgesWeight2or3 * 0.5;
    double uniquenessScore = analysis.uniqueNodesUsed - analysis.repeatNodeUsages * 0.5;
    
    return pathScore + uniquenessScore;
}

double Fitness::calculateFitness(
    const std::shared_ptr<std::vector<int>>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) const {
    if (!individual || individual->empty()) {
        return 0.0;
    }

    // Convert individual to DNA sequence
    std::vector<char> dna;
    for (int gene : *individual) {
        if (gene >= 0 && gene < static_cast<int>(instance.getSpectrum().size())) {
            const auto& kmer = instance.getSpectrum()[gene];
            dna.insert(dna.end(), kmer.begin(), kmer.end());
        } else {
            LOG_ERROR("Invalid gene value: " + std::to_string(gene));
            return 0.0;
        }
    }

    return calculateDNAFitness(dna, instance);
}

std::vector<std::vector<PreprocessedEdge>> OptimizedGraphBasedFitness::buildAdjacencyMatrix(
    const std::vector<std::vector<Edge>>& graph) const {
    
    const size_t n = graph.size();
    std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix(n, std::vector<PreprocessedEdge>(n));

    // Inicjalizacja macierzy sąsiedztwa
    for (size_t i = 0; i < n; ++i) {
        for (const auto& edge : graph[i]) {
            adjacencyMatrix[i][edge.to] = PreprocessedEdge(edge.to, edge.weight, true);
        }
    }

    return adjacencyMatrix;
}
