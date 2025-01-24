//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef FITNESS_H
#define FITNESS_H

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
#include <omp.h>

/**
 * Bardzo prosta funkcja fitness – liczy liczbę k-merów
 * wspólnych z oryginalnym spektrum DNA.
 */
class SimpleFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& individual,
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
        const std::shared_ptr<std::vector<int>>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;
};

/**
 * Przykład użycia algorytmu Smith-Waterman do oceny podobieństwa.
 */
class SmithWatermanFitness : public IFitness {
private:
    static constexpr int MATCH_SCORE = 5;
    static constexpr int MISMATCH_SCORE = -5;
    static constexpr int GAP_PENALTY = -5;
    static constexpr int GAP_EXTENSION = -5;

    int smithWaterman(const std::string& seq1, const std::string& seq2) const;

public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;
};

class OptimizedGraphBasedFitness : public IFitness {
public:
    struct Edge {
        int to;
        int weight;
        Edge(int t, int w) : to(t), weight(w) {}
    };
    
    struct PathAnalysis {
        int edgesWeight1;
        int edgesWeight2or3;
        int uniqueNodesUsed;
        int repeatNodeUsages;
        std::vector<int>& nodeUsageCount;
        
        PathAnalysis(std::vector<int>& buffer) 
            : edgesWeight1(0), edgesWeight2or3(0), 
              uniqueNodesUsed(0), repeatNodeUsages(0),
              nodeUsageCount(buffer) {}
    };

    // Główna funkcja budująca graf (z ewentualnym keszowaniem).
    std::vector<std::vector<Edge>> buildSpectrumGraph(
        const std::vector<std::string>& spectrum, 
        int k) const;
    
    // Analiza ścieżki w grafie – zliczamy wagi i liczbę użyć węzłów
    PathAnalysis analyzePath(const std::vector<int>& path,
                             const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix) const;
    void initBuffers(size_t size) const;
    OptimizedGraphBasedFitness() = default;
    
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;

    // Gdybyśmy chcieli wyczyścić kesz
    void clearCache() {
        graphCache.clear();
        nodeUsageBuffer.clear();
    }

private:
    // Bufor na liczbę użyć węzłów, by nie alokować za każdym razem
    mutable std::vector<int> nodeUsageBuffer;

    // Keszowanie wyników budowy grafu (by nie przeliczać 2x) 
    // Key: "k + ':' + posklejane_kmery", Value: tablica list sąsiedztwa
    mutable std::unordered_map<std::string, std::vector<std::vector<Edge>>> graphCache;

    std::string createCacheKey(const std::vector<std::string>& spectrum, int k) const;
    
    // Pomocnicza funkcja do ustalania wagi krawędzi
    int calculateEdgeWeight(
        const std::string& from, 
        const std::string& to, 
        int k) const;
    
    // Konwersja permutacji na "ścieżkę"
    const std::vector<int>& permutationToPath(std::shared_ptr<std::vector<int>> individual) const;

    // Budowa macierzy sąsiedztwa z grafu
    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(
        const std::vector<std::vector<Edge>>& graph) const;
};

class Fitness : public IFitness {
private:
    std::shared_ptr<IPopulationCache> m_cache;

protected:
    double calculateDNAFitness(const std::vector<char>& dna, const DNAInstance& instance) const;

public:
    explicit Fitness(std::shared_ptr<IPopulationCache> cache = nullptr) : m_cache(cache) {}
    
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override;
};

#endif //FITNESS_H
