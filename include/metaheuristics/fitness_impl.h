#pragma once

#include "../interfaces/i_fitness.h"
#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>
#include "preprocessed_edge.h"
#include "individual.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include <cmath>

class SimpleFitness : public IFitness {
public:
    SimpleFitness() : m_debugMode(false) {}
    
    void setDebugMode(bool debug) { m_debugMode = debug; }
    
    double calculateFitness(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

    double calculateConnectivityScore(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const;

    double calculateSpectrumCoverageScore(
        const std::vector<char>& dna,
        const DNAInstance& instance) const;

    double calculateLengthPenalty(
        int actualLength,
        int targetLength) const;

    int calculateLevenshteinDistance(
        const std::string& s1,
        const std::string& s2) const;

protected:
    virtual double calculateCoverage(const std::shared_ptr<Individual>& solution,
                                   const DNAInstance& instance) const = 0;

    virtual double calculateConnectivity(const std::shared_ptr<Individual>& solution,
                                       const DNAInstance& instance) const = 0;

private:
    bool m_debugMode;
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

class OptimizedGraphBasedFitness : public SimpleFitness {
private:
    GAConfig m_config;  // Changed from reference to value
    std::vector<std::shared_ptr<Individual>> m_currentPopulation;
    
    double calculateDistance(const std::shared_ptr<Individual>& ind1, 
                           const std::shared_ptr<Individual>& ind2) const {
        const auto& genes1 = ind1->getGenes();
        const auto& genes2 = ind2->getGenes();
        
        if (genes1.size() != genes2.size()) return 1.0;
        
        int differences = 0;
        for (size_t i = 0; i < genes1.size(); i++) {
            if (genes1[i] != genes2[i]) differences++;
        }
        
        return static_cast<double>(differences) / genes1.size();
    }
    
    double calculateSharingFactor(double distance, double radius, double alpha) const {
        if (distance >= radius) return 0.0;
        return 1.0 - std::pow(distance / radius, alpha);
    }

public:
    OptimizedGraphBasedFitness() : m_config(GAConfig()) {}  // Default constructor
    explicit OptimizedGraphBasedFitness(const GAConfig& config) 
        : m_config(config) {}

    void setCurrentPopulation(const std::vector<std::shared_ptr<Individual>>& population) {
        m_currentPopulation = population;
    }

    const std::vector<std::shared_ptr<Individual>>& getCurrentPopulation() const {
        return m_currentPopulation;
    }

    double calculateEdgeQuality(const std::shared_ptr<Individual>& individual,
                              const DNAInstance& instance) const;

    double calculateLength(const std::shared_ptr<Individual>& individual,
                         const DNAInstance& instance) const;

    double calculateFitness(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const override;

    double calculateConnectivityScore(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const;

    double calculateSpectrumCoverageScore(
        const std::vector<char>& dna,
        const DNAInstance& instance) const;

    double calculateLengthPenalty(
        int actualLength,
        int targetLength) const;

    int calculateLevenshteinDistance(
        const std::string& s1,
        const std::string& s2) const;

protected:
    double calculateCoverage(const std::shared_ptr<Individual>& solution,
                           const DNAInstance& instance) const override;
                           
    double calculateConnectivity(const std::shared_ptr<Individual>& solution,
                               const DNAInstance& instance) const override;

    std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(
        const DNAInstance& instance) const;
}; 