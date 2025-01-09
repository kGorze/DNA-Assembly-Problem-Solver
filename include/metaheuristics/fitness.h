//
// Created by konrad_guest on 07/01/2025.
//

#ifndef FITNESS_H
#define FITNESS_H

#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"
#include "utils/utility_functions.h"

#include <vector>
#include <unordered_map>

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

class GraphBasedFitness : public IFitness {
private:
    // Structure to represent graph edges
    struct Edge {
        int to;
        int weight;
        Edge(int t, int w) : to(t), weight(w) {}
    };
    
    // Helper struct to store path analysis results
    struct PathAnalysis {
        int edgesWeight1 = 0;
        int edgesWeight2or3 = 0;
        int uniqueNodesUsed = 0;
        int repeatNodeUsages = 0;
        std::unordered_map<int, int> nodeUsageCount;
    };

    // Helper methods
    std::vector<std::vector<Edge>> buildSpectrumGraph(const std::vector<std::string>& spectrum, int k) const;
    int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;
    PathAnalysis analyzePath(const std::vector<int>& path, 
                           const std::vector<std::vector<Edge>>& graph) const;
    std::vector<int> permutationToPath(void* individual) const;

public:
    double evaluate(void* individual,
                   const DNAInstance& instance,
                   std::shared_ptr<IRepresentation> representation) const override;
};

#endif //FITNESS_H
