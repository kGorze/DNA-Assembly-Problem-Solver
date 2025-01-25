#pragma once

#include "../interfaces/i_crossover.h"
#include "../dna/dna_instance.h"
#include "individual.h"
#include <vector>
#include <memory>
#include <random>
#include <algorithm>

class OnePointCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (parents.size() < 2) {
            throw std::invalid_argument("Need at least 2 parents for crossover");
        }

        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);

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

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(1, static_cast<int>(genes1.size() - 1));
        int crossPoint = dis(gen);

        // Create first offspring
        std::vector<int> offspring1Genes(genes1.begin(), genes1.begin() + crossPoint);
        offspring1Genes.insert(offspring1Genes.end(), genes2.begin() + crossPoint, genes2.end());
        offspring.push_back(std::make_shared<Individual>(offspring1Genes));

        // Create second offspring
        std::vector<int> offspring2Genes(genes2.begin(), genes2.begin() + crossPoint);
        offspring2Genes.insert(offspring2Genes.end(), genes1.begin() + crossPoint, genes1.end());
        offspring.push_back(std::make_shared<Individual>(offspring2Genes));

        return offspring;
    }
};

class OrderCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
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

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, static_cast<int>(genes1.size() - 1));
        
        int start = dis(gen);
        int end = dis(gen);
        if (start > end) std::swap(start, end);

        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);

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
        for (int i = 0; i < genes1.size(); i++) {
            if (i >= start && i <= end) continue;
            while (j < genes2.size() && used1[genes2[j]]) j++;
            if (j < genes2.size()) {
                offspring1Genes[i] = genes2[j++];
            }
        }
        
        offspring.push_back(std::make_shared<Individual>(offspring1Genes));

        // Create second offspring (reverse roles)
        std::vector<int> offspring2Genes(genes1.size());
        std::vector<bool> used2(genes1.size(), false);
        
        for (int i = start; i <= end; i++) {
            offspring2Genes[i] = genes2[i];
            used2[genes2[i]] = true;
        }
        
        j = 0;
        for (int i = 0; i < genes1.size(); i++) {
            if (i >= start && i <= end) continue;
            while (j < genes1.size() && used2[genes1[j]]) j++;
            if (j < genes1.size()) {
                offspring2Genes[i] = genes1[j++];
            }
        }
        
        offspring.push_back(std::make_shared<Individual>(offspring2Genes));

        return offspring;
    }
};

class EdgeRecombination : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        
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

        std::random_device rd;
        std::mt19937 gen(rd());

        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);

        // Create two offspring
        for (int k = 0; k < 2; k++) {
            std::vector<int> offspringGenes;
            offspringGenes.reserve(genes1.size());
            std::vector<bool> used(genes1.size(), false);

            // Start with a random gene
            std::uniform_int_distribution<> dis(0, static_cast<int>(genes1.size() - 1));
            int current = dis(gen);
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
                        if (!used[i]) unused.push_back(i);
                    }
                    if (!unused.empty()) {
                        std::uniform_int_distribution<> dis2(0, static_cast<int>(unused.size() - 1));
                        next = unused[dis2(gen)];
                    }
                }
                
                if (next != -1) {
                    offspringGenes.push_back(next);
                    used[next] = true;
                    current = next;
                }
            }
            
            offspring.push_back(std::make_shared<Individual>(offspringGenes));
        }

        return offspring;
    }
};

class PMXCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
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
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
}; 