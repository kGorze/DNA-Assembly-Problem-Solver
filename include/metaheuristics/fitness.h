//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef FITNESS_H
#define FITNESS_H

#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"
#include "utils/utility_functions.h"

#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <omp.h>

/**
 * Interfejs fitness
 */
class IFitness {
public:
    virtual ~IFitness() = default;
    /**
     * Ocena przystosowania pojedynczego osobnika
     */
    virtual double evaluate(std::shared_ptr<std::vector<int>> individual,
                            const DNAInstance &instance,
                            std::shared_ptr<IRepresentation> representation) const = 0;
};

/**
 * Bardzo prosta funkcja fitness – liczy liczbę k-merów
 * wspólnych z oryginalnym spektrum DNA.
 */
class SimpleFitness : public IFitness {
public:
    double evaluate(std::shared_ptr<std::vector<int>> individual,
                    const DNAInstance &instance,
                    std::shared_ptr<IRepresentation> representation) const override;
};

/**
 * Ulepszona wersja fitness – oprócz liczby trafionych k-merów
 * uwzględnia odległość Levenshteina do oryginalnego DNA.
 */
class BetterFitness : public IFitness {
public:
    double evaluate(std::shared_ptr<std::vector<int>> individual,
                    const DNAInstance &instance,
                    std::shared_ptr<IRepresentation> representation) const override;
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
    double evaluate(std::shared_ptr<std::vector<int>> individual, 
                    const DNAInstance& instance,
                    std::shared_ptr<IRepresentation> representation) const override;
};

/**
 * Bardziej zaawansowana wersja fitness, wykorzystująca "graf" k-merów
 * i analizę ścieżki w tym grafie.
 *
 * Poprawiona wersja, w której macierz adjacencyMatrix tworzona jest
 * lokalnie w metodzie evaluate(...) aby uniknąć błędów wątkowych.
 */
class OptimizedGraphBasedFitness : public IFitness {
private:
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

    // Bufor na liczbę użyć węzłów, by nie alokować za każdym razem
    mutable std::vector<int> nodeUsageBuffer;

    // Keszowanie wyników budowy grafu (by nie przeliczać 2x) 
    // Key: "k + ':' + posklejane_kmery", Value: tablica list sąsiedztwa
    mutable std::unordered_map<std::string, std::vector<std::vector<Edge>>> graphCache;

    // Usuwamy 'mutable std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix;'
    // bo powoduje problem we współbieżności i mamy teraz lokalną alokację w evaluate().

    
    std::string createCacheKey(const std::vector<std::string>& spectrum, int k) const;

    
    // Pomocnicza funkcja do ustalania wagi krawędzi
    int calculateEdgeWeight(
        const std::string& from, 
        const std::string& to, 
        int k) const;
    
 
    
    // Konwersja permutacji na "ścieżkę"
    const std::vector<int>& permutationToPath(std::shared_ptr<std::vector<int>> individual) const;

public:

        
    // Główna funkcja budująca graf (z ewentualnym keszowaniem).
    std::vector<std::vector<Edge>> buildSpectrumGraph(
        const std::vector<std::string>& spectrum, 
        int k) const;
    
    struct PreprocessedEdge {
        int to;
        int weight;
        bool exists;
        PreprocessedEdge(int t = 0, int w = 0, bool e = false) 
            : to(t), weight(w), exists(e) {}
    };
    
    // Analiza ścieżki w grafie – zliczamy wagi i liczbę użyć węzłów
    PathAnalysis analyzePath(const std::vector<int>& path,
                             const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix) const;
    void initBuffers(size_t size) const;
    OptimizedGraphBasedFitness() = default;
    
    double evaluate(std::shared_ptr<std::vector<int>> individual,
                    const DNAInstance& instance,
                    std::shared_ptr<IRepresentation> representation) const override;

    // Gdybyśmy chcieli wyczyścić kesz
    void clearCache() {
        graphCache.clear();
        nodeUsageBuffer.clear();
    }
};

#endif //FITNESS_H
