//
// Created by konrad_guest on 07/01/2025.
//

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

class IFitness {
public:
    virtual ~IFitness() = default;
    /**
     * Ocena przystosowania pojedynczego osobnika
     */
    virtual double evaluate(void* individual,
                        const DNAInstance &instance,
                        std::shared_ptr<IRepresentation> representation) const = 0;
};

/**
 * Przykładowa implementacja – do rozbudowy
 */
class SimpleFitness : public IFitness {
public:
    double evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const override;
};

class BetterFitness : public IFitness {
public:
    double evaluate(void* individual,
                    const DNAInstance &instance,
                    std::shared_ptr<IRepresentation> representation) const override;
};

class SmithWatermanFitness : public IFitness {
private:
    // Smith-Waterman scoring parameters
    static constexpr int MATCH_SCORE = 5;      // Will be divided by 5 later
    static constexpr int MISMATCH_SCORE = -5;  // Will be divided by 5 later
    static constexpr int GAP_PENALTY = -5;     // Gap insertion and extension
    static constexpr int GAP_EXTENSION = -5;   // Same as insertion in this case

    int smithWaterman(const std::string& seq1, const std::string& seq2) const;

public:
    double evaluate(void* individual, 
                   const DNAInstance& instance,
                   std::shared_ptr<IRepresentation> representation) const override;
};

class OptimizedGraphBasedFitness : public IFitness {
private:
    // Struktury pomocnicze
    struct Edge {
        int to;
        int weight;
        Edge(int t, int w) : to(t), weight(w) {}
    };
    
    struct PreprocessedEdge {
        int to;
        int weight;
        bool exists;
        PreprocessedEdge(int t = 0, int w = 0, bool e = false) 
            : to(t), weight(w), exists(e) {}
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

    // Bufory i cache
    mutable std::vector<int> nodeUsageBuffer;
    mutable std::unordered_map<std::string, std::vector<std::vector<Edge>>> graphCache;
    mutable std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix;

    // Metody pomocnicze
    void initBuffers(size_t size) const;
    void preprocessGraph(const std::vector<std::vector<Edge>>& graph) const;
    std::string createCacheKey(const std::vector<std::string>& spectrum, int k) const;
    
    std::vector<std::vector<Edge>> buildSpectrumGraph(
        const std::vector<std::string>& spectrum, 
        int k) const;
    
    int calculateEdgeWeight(
        const std::string& from, 
        const std::string& to, 
        int k) const;
    
    PathAnalysis analyzePath(
        const std::vector<int>& path,
        size_t graphSize) const;
    
    const std::vector<int>& permutationToPath(void* individual) const;

public:
    OptimizedGraphBasedFitness() = default;
    
    double evaluate(
        void* individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

    // Czyszczenie cache
    void clearCache() {
        graphCache.clear();
        adjacencyMatrix.clear();
        nodeUsageBuffer.clear();
    }
};


#endif //FITNESS_H
