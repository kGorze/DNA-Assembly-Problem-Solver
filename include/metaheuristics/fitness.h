//
// Created by konrad_guest on 07/01/2025.
// SMART

#pragma once

#include "../interfaces/i_representation.h"
#include "../interfaces/i_fitness.h"
#include "../interfaces/i_population_cache.h"
#include "../generator/dna_generator.h"
#include "../utils/utility_functions.h"
#include "../metaheuristics/path_analyzer.h"
#include "../metaheuristics/preprocessed_edge.h"
#include "../metaheuristics/individual.h"
#include "../utils/logger.h"
#include "../dna/dna_instance.h"

#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include <mutex>
#include <cmath>
#include <limits>

/**
 * Simple fitness function that counts the number of k-mers
 * shared with the original DNA spectrum.
 */
class SimpleFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override {
        if (!individual || !representation) return -std::numeric_limits<double>::infinity();
        
        try {
            const auto& genes = individual->getGenes();
            const auto& spectrum = instance.getSpectrum();
            
            if (genes.empty()) return 0.0;
            
            int matches = 0;
            for (size_t i = 0; i < genes.size(); ++i) {
                if (genes[i] >= 0 && genes[i] < static_cast<int>(spectrum.size())) {
                    matches++;
                }
            }
            
            return static_cast<double>(matches);
        } catch (const std::exception& e) {
            LOG_ERROR("Error calculating simple fitness: " + std::string(e.what()));
            return -std::numeric_limits<double>::infinity();
        }
    }
};

/**
 * Enhanced fitness function that considers both k-mer matches
 * and Levenshtein distance to the original DNA.
 */
class BetterFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override {
        if (!individual || !representation) return -std::numeric_limits<double>::infinity();
        
        try {
            const auto& genes = individual->getGenes();
            const auto& spectrum = instance.getSpectrum();
            
            if (genes.empty()) return 0.0;
            
            // Calculate k-mer matches
            int matches = 0;
            for (size_t i = 0; i < genes.size(); ++i) {
                if (genes[i] >= 0 && genes[i] < static_cast<int>(spectrum.size())) {
                    matches++;
                }
            }
            
            // Calculate Levenshtein distance
            const std::string dna = representation->toString(genes);
            const std::string original = instance.getOriginalDNA();
            const double distance = levenshteinDistance(dna, original);
            
            // Combine scores
            const double matchScore = static_cast<double>(matches);
            const double distanceScore = 1.0 / (1.0 + distance);
            
            return matchScore * distanceScore;
        } catch (const std::exception& e) {
            LOG_ERROR("Error calculating better fitness: " + std::string(e.what()));
            return -std::numeric_limits<double>::infinity();
        }
    }

private:
    double levenshteinDistance(const std::string& s1, const std::string& s2) const {
        if (s1.empty()) return s2.length();
        if (s2.empty()) return s1.length();
        
        std::vector<std::vector<double>> dp(s1.length() + 1, 
                                          std::vector<double>(s2.length() + 1));
        
        for (size_t i = 0; i <= s1.length(); i++) dp[i][0] = i;
        for (size_t j = 0; j <= s2.length(); j++) dp[0][j] = j;
        
        for (size_t i = 1; i <= s1.length(); i++) {
            for (size_t j = 1; j <= s2.length(); j++) {
                dp[i][j] = std::min({
                    dp[i-1][j] + 1,
                    dp[i][j-1] + 1,
                    dp[i-1][j-1] + (s1[i-1] != s2[j-1])
                });
            }
        }
        
        return dp[s1.length()][s2.length()];
    }
};

/**
 * Smith-Waterman algorithm based fitness function.
 */
class SmithWatermanFitness : public IFitness {
public:
    static constexpr int MATCH_SCORE = 2;
    static constexpr int MISMATCH_SCORE = -1;
    static constexpr int GAP_PENALTY = -2;

    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override {
        if (!individual || !representation) return -std::numeric_limits<double>::infinity();
        
        try {
            const auto& genes = individual->getGenes();
            if (genes.empty()) return 0.0;
            
            const std::string dna = representation->toString(genes);
            const std::string original = instance.getOriginalDNA();
            
            return static_cast<double>(smithWaterman(dna, original));
        } catch (const std::exception& e) {
            LOG_ERROR("Error calculating Smith-Waterman fitness: " + std::string(e.what()));
            return -std::numeric_limits<double>::infinity();
        }
    }

private:
    int smithWaterman(const std::string& seq1, const std::string& seq2) const {
        if (seq1.empty() || seq2.empty()) return 0;
        
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
};

/**
 * Graph-based fitness function with optimized caching.
 */
class OptimizedGraphBasedFitness : public IFitness {
public:
    OptimizedGraphBasedFitness() = default;
    ~OptimizedGraphBasedFitness() override = default;
    
    double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const override {
        if (!individual || !representation) return -std::numeric_limits<double>::infinity();
        
        try {
            const auto& genes = individual->getGenes();
            if (genes.empty()) return 0.0;
            
            const std::string dna = representation->toString(genes);
            
            // Check cache first
            {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                auto it = m_cache.find(dna);
                if (it != m_cache.end()) {
                    return it->second;
                }
            }
            
            // Calculate new fitness
            double fitness = 0.0;
            {
                std::lock_guard<std::mutex> lock(m_bufferMutex);
                ensureBufferSizes(genes.size());
                buildGraph(genes, instance);
                
                const double connectivityScore = calculateConnectivityScore();
                const double coverageScore = calculateSpectrumCoverageScore(genes, instance);
                const double lengthPenalty = calculateLengthPenalty(genes, instance);
                
                fitness = connectivityScore * coverageScore * (1.0 - lengthPenalty);
            }
            
            // Cache the result
            {
                std::lock_guard<std::mutex> lock(m_cacheMutex);
                if (m_cache.size() >= maxCacheSize) {
                    cleanupCache();
                }
                m_cache[dna] = fitness;
            }
            
            return fitness;
        } catch (const std::exception& e) {
            LOG_ERROR("Error calculating graph-based fitness: " + std::string(e.what()));
            return -std::numeric_limits<double>::infinity();
        }
    }

private:
    static constexpr size_t maxCacheSize = 10000;
    
    mutable std::mutex m_cacheMutex;
    mutable std::unordered_map<std::string, double> m_cache;
    
    mutable std::mutex m_bufferMutex;
    mutable std::vector<int> m_nodeUsageBuffer;
    mutable std::vector<bool> m_coveredBuffer;
    mutable std::vector<std::vector<PreprocessedEdge>> m_adjacencyMatrixBuffer;
    
    void ensureBufferSizes(size_t size) const {
        m_nodeUsageBuffer.resize(size);
        m_coveredBuffer.resize(size);
        m_adjacencyMatrixBuffer.resize(size, std::vector<PreprocessedEdge>(size));
    }
    
    void buildGraph(const std::vector<int>& genes, const DNAInstance& instance) const {
        // Implementation details...
    }
    
    bool areNodesConnected(int node1, int node2, const DNAInstance& instance) const {
        // Implementation details...
        return false;
    }
    
    double calculateConnectivityScore() const {
        // Implementation details...
        return 0.0;
    }
    
    double calculateSpectrumCoverageScore(
        const std::vector<int>& genes, 
        const DNAInstance& instance) const {
        // Implementation details...
        return 0.0;
    }
    
    double calculateLengthPenalty(
        const std::vector<int>& genes,
        const DNAInstance& instance) const {
        // Implementation details...
        return 0.0;
    }
    
    void cleanupCache() const {
        // Remove oldest entries when cache is full
        if (m_cache.size() > maxCacheSize / 2) {
            m_cache.clear();
        }
    }
};
