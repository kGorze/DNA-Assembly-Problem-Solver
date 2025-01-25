#pragma once

#include "../interfaces/i_fitness.h"
#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>

class SimpleFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

private:
    double calculateConnectivityScore(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const;

    double calculateSpectrumCoverageScore(
        const std::vector<char>& dna,
        const DNAInstance& instance) const;

    double calculateLengthPenalty(
        int actualLength,
        int targetLength) const;
};

class BetterFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

private:
    double calculateConnectivityScore(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const;

    double calculateSpectrumCoverageScore(
        const std::vector<char>& dna,
        const DNAInstance& instance) const;

    double calculateLengthPenalty(
        int actualLength,
        int targetLength) const;

    double calculateSequenceSimilarity(
        const std::vector<char>& dna,
        const DNAInstance& instance) const;
};

class SmithWatermanFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

private:
    double smithWaterman(const std::string& seq1, const std::string& seq2) const;
};

struct PreprocessedEdge {
    int to;
    double weight;
    bool exists;

    PreprocessedEdge() : to(-1), weight(0.0), exists(false) {}
    PreprocessedEdge(int t, double w, bool e) : to(t), weight(w), exists(e) {}
};

class OptimizedGraphBasedFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

private:
    double calculateConnectivityScore(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const;

    double calculateSpectrumCoverageScore(
        const std::vector<char>& dna,
        const DNAInstance& instance) const;

    double calculateLengthPenalty(
        int actualLength,
        int targetLength) const;

    double calculatePathOptimalityScore(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const;

    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(
        const DNAInstance& instance) const;

    int calculateEdgeWeight(
        const std::string& from,
        const std::string& to,
        int k) const;
}; 