#pragma once

#include "../interfaces/i_crossover.h"
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
    std::pair<std::vector<int>, std::vector<int>> validateAndGetGenes(
        const std::vector<std::shared_ptr<Individual>>& parents) {
        
        if (parents.size() < 2) {
            throw std::invalid_argument("Need at least 2 parents for crossover");
        }

        const auto& parent1 = parents[0];
        const auto& parent2 = parents[1];

        if (!parent1 || !parent2) {
            throw std::invalid_argument("Invalid parent pointers");
        }

        const auto& genes1 = parent1->getGenes();
        const auto& genes2 = parent2->getGenes();

        if (genes1.size() != genes2.size()) {
            throw std::invalid_argument("Parents must have same length");
        }
        
        return {genes1, genes2};
    }
    
    // Create offspring with validation
    std::shared_ptr<Individual> createOffspring(
        std::vector<int> genes,
        std::shared_ptr<IRepresentation> representation,
        const DNAInstance& instance) {
        
        try {
            auto offspring = std::make_shared<Individual>(std::move(genes));
            if (!offspring) {
                LOG_ERROR("Failed to create offspring");
                return nullptr;
            }
            
            if (!representation->isValid(offspring, instance)) {
                LOG_WARNING("Created invalid offspring");
                return nullptr;
            }
            
            return offspring;
        } catch (const std::exception& e) {
            LOG_ERROR("Error creating offspring: " + std::string(e.what()));
            return nullptr;
        }
    }
}

class OnePointCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (!representation) {
            throw std::invalid_argument("Invalid representation pointer");
        }
        
        try {
            auto [genes1, genes2] = validateAndGetGenes(parents);
            
            std::vector<std::shared_ptr<Individual>> offspring;
            offspring.reserve(2);
            
            auto& rng = Random::instance();
            int crossPoint = rng.getRandomInt(1, static_cast<int>(genes1.size() - 1));
            
            // Create first offspring
            std::vector<int> offspring1Genes(genes1.begin(), genes1.begin() + crossPoint);
            offspring1Genes.insert(offspring1Genes.end(), genes2.begin() + crossPoint, genes2.end());
            
            if (auto child = createOffspring(offspring1Genes, representation, instance)) {
                offspring.push_back(child);
            }
            
            // Create second offspring
            std::vector<int> offspring2Genes(genes2.begin(), genes2.begin() + crossPoint);
            offspring2Genes.insert(offspring2Genes.end(), genes1.begin() + crossPoint, genes1.end());
            
            if (auto child = createOffspring(offspring2Genes, representation, instance)) {
                offspring.push_back(child);
            }
            
            return offspring;
            
        } catch (const std::exception& e) {
            LOG_ERROR("Error in one-point crossover: " + std::string(e.what()));
            return {};
        }
    }
};

class OrderCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (!representation) {
            throw std::invalid_argument("Invalid representation pointer");
        }
        
        try {
            auto [genes1, genes2] = validateAndGetGenes(parents);
            
            std::vector<std::shared_ptr<Individual>> offspring;
            offspring.reserve(2);
            
            auto& rng = Random::instance();
            int start = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
            int end = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
            if (start > end) std::swap(start, end);
            
            // Create first offspring
            std::vector<int> offspring1Genes(genes1.size());
            std::vector<bool> used1(genes1.size(), false);
            
            // Copy segment from first parent
            for (int i = start; i <= end; i++) {
                offspring1Genes[i] = genes1[i];
                used1[genes1[i]] = true;
            }
            
            // Fill remaining positions from second parent
            int j = 0;
            for (int i = 0; i < static_cast<int>(genes1.size()); i++) {
                if (i >= start && i <= end) continue;
                while (j < static_cast<int>(genes2.size()) && used1[genes2[j]]) j++;
                if (j < static_cast<int>(genes2.size())) {
                    offspring1Genes[i] = genes2[j++];
                }
            }
            
            if (auto child = createOffspring(offspring1Genes, representation, instance)) {
                offspring.push_back(child);
            }
            
            // Create second offspring (reverse roles)
            std::vector<int> offspring2Genes(genes1.size());
            std::vector<bool> used2(genes1.size(), false);
            
            for (int i = start; i <= end; i++) {
                offspring2Genes[i] = genes2[i];
                used2[genes2[i]] = true;
            }
            
            j = 0;
            for (int i = 0; i < static_cast<int>(genes1.size()); i++) {
                if (i >= start && i <= end) continue;
                while (j < static_cast<int>(genes1.size()) && used2[genes1[j]]) j++;
                if (j < static_cast<int>(genes1.size())) {
                    offspring2Genes[i] = genes1[j++];
                }
            }
            
            if (auto child = createOffspring(offspring2Genes, representation, instance)) {
                offspring.push_back(child);
            }
            
            return offspring;
            
        } catch (const std::exception& e) {
            LOG_ERROR("Error in order crossover: " + std::string(e.what()));
            return {};
        }
    }
};

class EdgeRecombination : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (!representation) {
            throw std::invalid_argument("Invalid representation pointer");
        }
        
        try {
            auto [genes1, genes2] = validateAndGetGenes(parents);
            
            std::vector<std::shared_ptr<Individual>> offspring;
            offspring.reserve(2);
            
            // Build edge table
            std::vector<std::vector<int>> edgeTable(genes1.size());
            for (size_t i = 0; i < genes1.size(); i++) {
                // Add edges from first parent
                edgeTable[genes1[i]].push_back(genes1[(i + 1) % genes1.size()]);
                edgeTable[genes1[i]].push_back(genes1[(i + genes1.size() - 1) % genes1.size()]);
                
                // Add edges from second parent
                edgeTable[genes2[i]].push_back(genes2[(i + 1) % genes2.size()]);
                edgeTable[genes2[i]].push_back(genes2[(i + genes2.size() - 1) % genes2.size()]);
                
                // Remove duplicates
                std::sort(edgeTable[genes1[i]].begin(), edgeTable[genes1[i]].end());
                edgeTable[genes1[i]].erase(
                    std::unique(edgeTable[genes1[i]].begin(), edgeTable[genes1[i]].end()),
                    edgeTable[genes1[i]].end()
                );
            }
            
            auto& rng = Random::instance();
            
            // Create two offspring
            for (int k = 0; k < 2; k++) {
                std::vector<int> offspringGenes;
                offspringGenes.reserve(genes1.size());
                std::vector<bool> used(genes1.size(), false);
                
                // Start with a random gene
                int current = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
                offspringGenes.push_back(current);
                used[current] = true;
                
                // Build the rest of the path
                while (offspringGenes.size() < genes1.size()) {
                    const auto& neighbors = edgeTable[current];
                    
                    // Find unused neighbor with fewest remaining neighbors
                    int next = -1;
                    int minNeighbors = std::numeric_limits<int>::max();
                    
                    for (int neighbor : neighbors) {
                        if (used[neighbor]) continue;
                        
                        int remainingNeighbors = 0;
                        for (int n : edgeTable[neighbor]) {
                            if (!used[n]) remainingNeighbors++;
                        }
                        
                        if (remainingNeighbors < minNeighbors) {
                            minNeighbors = remainingNeighbors;
                            next = neighbor;
                        }
                    }
                    
                    // If no unused neighbors, pick random unused gene
                    if (next == -1) {
                        std::vector<int> unused;
                        for (size_t i = 0; i < used.size(); i++) {
                            if (!used[i]) unused.push_back(static_cast<int>(i));
                        }
                        if (!unused.empty()) {
                            next = unused[rng.getRandomInt(0, static_cast<int>(unused.size() - 1))];
                        }
                    }
                    
                    if (next != -1) {
                        offspringGenes.push_back(next);
                        used[next] = true;
                        current = next;
                    }
                }
                
                if (auto child = createOffspring(offspringGenes, representation, instance)) {
                    offspring.push_back(child);
                }
            }
            
            return offspring;
            
        } catch (const std::exception& e) {
            LOG_ERROR("Error in edge recombination crossover: " + std::string(e.what()));
            return {};
        }
    }
};

class PMXCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (!representation) {
            throw std::invalid_argument("Invalid representation pointer");
        }
        
        try {
            auto [genes1, genes2] = validateAndGetGenes(parents);
            
            std::vector<std::shared_ptr<Individual>> offspring;
            offspring.reserve(2);
            
            auto& rng = Random::instance();
            int start = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
            int end = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
            if (start > end) std::swap(start, end);
            
            // Create first offspring
            std::vector<int> offspring1Genes(genes1.size(), -1);
            std::vector<bool> used1(genes1.size(), false);
            
            // Copy segment from first parent
            for (int i = start; i <= end; i++) {
                offspring1Genes[i] = genes1[i];
                used1[genes1[i]] = true;
            }
            
            // Map corresponding segment from second parent
            std::unordered_map<int, int> mapping;
            for (int i = start; i <= end; i++) {
                if (genes1[i] != genes2[i]) {
                    mapping[genes2[i]] = genes1[i];
                }
            }
            
            // Fill remaining positions
            for (int i = 0; i < static_cast<int>(genes1.size()); i++) {
                if (i >= start && i <= end) continue;
                
                int value = genes2[i];
                while (mapping.find(value) != mapping.end()) {
                    value = mapping[value];
                }
                offspring1Genes[i] = value;
            }
            
            if (auto child = createOffspring(offspring1Genes, representation, instance)) {
                offspring.push_back(child);
            }
            
            // Create second offspring (reverse roles)
            std::vector<int> offspring2Genes(genes1.size(), -1);
            std::vector<bool> used2(genes1.size(), false);
            
            for (int i = start; i <= end; i++) {
                offspring2Genes[i] = genes2[i];
                used2[genes2[i]] = true;
            }
            
            mapping.clear();
            for (int i = start; i <= end; i++) {
                if (genes1[i] != genes2[i]) {
                    mapping[genes1[i]] = genes2[i];
                }
            }
            
            for (int i = 0; i < static_cast<int>(genes1.size()); i++) {
                if (i >= start && i <= end) continue;
                
                int value = genes1[i];
                while (mapping.find(value) != mapping.end()) {
                    value = mapping[value];
                }
                offspring2Genes[i] = value;
            }
            
            if (auto child = createOffspring(offspring2Genes, representation, instance)) {
                offspring.push_back(child);
            }
            
            return offspring;
            
        } catch (const std::exception& e) {
            LOG_ERROR("Error in PMX crossover: " + std::string(e.what()));
            return {};
        }
    }
};

class DistancePreservingCrossover : public ICrossover {
private:
    class DistanceMatrix {
    public:
        explicit DistanceMatrix(const std::vector<int>& perm);
        int getDistance(int from, int to) const;
    private:
        std::vector<int> distances;
    };

public:
    std::vector<std::shared_ptr<Individual>> crossover(
        [[maybe_unused]] const std::vector<std::shared_ptr<Individual>>& parents,
        [[maybe_unused]] const DNAInstance& instance,
        [[maybe_unused]] std::shared_ptr<IRepresentation> representation) override {
        // Implementation here
        return std::vector<std::shared_ptr<Individual>>();
    }
}; 