#pragma once

#include "../interfaces/i_crossover.h"
#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include "individual.h"
#include "../utils/logging.h"
#include "utils/random.h"
#include "metaheuristics/dna_utils.h"
#include <vector>
#include <memory>
#include <random>
#include <algorithm>
#include <mutex>
#include <stdexcept>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
#include <queue>

namespace {
    // Validate parents and get their genes
    [[maybe_unused]]
    std::pair<std::vector<int>, std::vector<int>> validateAndGetGenes(
        const std::vector<std::shared_ptr<Individual>>& parents) {
        
        if (parents.size() < 2 || !parents[0] || !parents[1]) {
            return {};
        }

        const auto& parent1 = parents[0];
        const auto& parent2 = parents[1];

        const auto& genes1 = parent1->getGenes();
        const auto& genes2 = parent2->getGenes();

        if (genes1.size() != genes2.size()) {
            throw std::invalid_argument("Parents must have same length");
        }
        
        return {genes1, genes2};
    }
    
    // Create offspring with validation
    [[maybe_unused]]
    std::shared_ptr<Individual> createOffspring(
        std::vector<int> genes,
        std::shared_ptr<IRepresentation> representation,
        const DNAInstance& instance) {
        
        // Basic validation
        if (genes.empty()) {
            LOG_ERROR("Cannot create offspring with empty genes");
            return nullptr;
        }
        
        // Check if all genes are valid indices
        for (int gene : genes) {
            if (gene < 0 || static_cast<size_t>(gene) >= instance.getSpectrum().size()) {
                LOG_ERROR("Invalid gene index in offspring: " + std::to_string(gene));
                return nullptr;
            }
        }
        
        // Create individual and validate with representation
        auto offspring = std::make_shared<Individual>(std::move(genes));
        if (!representation->isValid(offspring, instance)) {
            LOG_ERROR("Created offspring failed representation validation");
            return nullptr;
        }
        
        return offspring;
    }
}

class OrderCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
};

class EdgeRecombination : public ICrossover {
public:
    struct EdgeInfo {
        int node;
        int overlapQuality;
        
        EdgeInfo(int n, int q) : node(n), overlapQuality(q) {}
        bool operator==(const EdgeInfo& other) const {
            return node == other.node;
        }
    };

    class EdgeTable {
    private:
        std::unordered_map<int, std::vector<EdgeInfo>> edges;

    public:
        EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents, const DNAInstance& instance);
        std::vector<EdgeInfo> getNeighbors(int node) const;  // Changed return type to EdgeInfo
        void removeNode(int node);
        bool hasNode(int node) const;
    };

    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
};

class PMXCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
};

class DNAAlignmentCrossover : public ICrossover {
private:
    struct AlignmentSegment {
        std::vector<int> genes;
        double quality;
        
        AlignmentSegment(std::vector<int> g, double q) : genes(std::move(g)), quality(q) {}
        
        bool operator<(const AlignmentSegment& other) const {
            return quality > other.quality;  // Sort by descending quality
        }
    };
    
    std::vector<AlignmentSegment> findAlignmentSegments(
        const std::vector<int>& genes,
        const DNAInstance& instance,
        int minLength = 3) const;
        
    double calculateSegmentQuality(
        const std::vector<int>& genes,
        size_t start,
        size_t length,
        const DNAInstance& instance) const;
        
    std::vector<int> mergeSegments(
        const std::vector<AlignmentSegment>& segments1,
        const std::vector<AlignmentSegment>& segments2,
        size_t targetLength) const;

public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override;
}; 