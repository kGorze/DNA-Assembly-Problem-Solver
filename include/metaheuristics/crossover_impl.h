#pragma once

#include "../interfaces/i_crossover.h"
#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include "individual.h"
#include "../utils/logging.h"
#include "utils/random.h"
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
        
        // Always create and return offspring, let fitness function handle validation
        return std::make_shared<Individual>(std::move(genes));
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
    class EdgeTable {
    private:
        std::unordered_map<int, std::vector<int>> edges;

    public:
        EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents);
        std::vector<int> getNeighbors(int node) const;
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