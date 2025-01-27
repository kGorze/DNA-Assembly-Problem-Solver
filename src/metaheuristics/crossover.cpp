//
// Created by konrad_guest on 07/01/2025.
// SMART

#include "../include/metaheuristics/crossover_impl.h"
#include "../include/utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <unordered_set>

/**
 * Helper function to verify if a chromosome is a valid permutation
 * in the range 0..(size-1).
 */
[[maybe_unused]]
static bool validateAndClampChromosome(std::shared_ptr<std::vector<int>> indiv, int size) {
    if (!indiv) return false;
    if ((int)indiv->size() != size) {
        return false;
    }
    std::vector<bool> used(size, false);
    for (int gene : *indiv) {
        if (gene < 0 || gene >= size || used[gene]) {
            return false;
        }
        used[gene] = true;
    }
    return true;
}

bool isValidPermutation(const std::vector<int>& perm, int size) {
    if (perm.size() != static_cast<size_t>(size)) return false;
    
    std::vector<bool> used(size, false);
    for (int val : perm) {
        if (val < 0 || val >= size || used[val]) {
            return false;
        }
        used[val] = true;
    }
    return true;
}

// ========================================
// =           ORDER CROSSOVER (OX)      =
// ========================================
std::vector<std::shared_ptr<Individual>> 
OrderCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (!representation) {
        LOG_ERROR("Invalid representation pointer");
        return parents;  // Return parents instead of empty vector
    }
    
    try {
        auto [genes1, genes2] = validateAndGetGenes(parents);
        if (genes1.empty() || genes2.empty()) {
            LOG_ERROR("Empty parent genes");
            return parents;  // Return parents instead of empty vector
        }
        
        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);
        
        auto& rng = Random::instance();
        int start = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
        int end = rng.getRandomInt(0, static_cast<int>(genes1.size() - 1));
        if (start > end) std::swap(start, end);
        
        int maxAttempts = 3;  // Try up to 3 times to create valid offspring
        while (offspring.size() < 2 && maxAttempts > 0) {
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
            bool validOffspring1 = true;
            for (int i = 0; i < static_cast<int>(genes1.size()); i++) {
                if (i >= start && i <= end) continue;
                
                // Find next unused gene from parent2
                while (j < static_cast<int>(genes2.size()) && used1[genes2[j]]) j++;
                
                if (j < static_cast<int>(genes2.size())) {
                    offspring1Genes[i] = genes2[j++];
                } else {
                    // If we run out of genes, find any unused valid index
                    bool foundUnused = false;
                    for (size_t val = 0; val < instance.getSpectrum().size(); val++) {
                        if (!used1[val]) {
                            offspring1Genes[i] = val;
                            used1[val] = true;
                            foundUnused = true;
                            break;
                        }
                    }
                    if (!foundUnused) {
                        validOffspring1 = false;
                        break;
                    }
                }
            }
            
            if (validOffspring1) {
                if (auto child = createOffspring(offspring1Genes, representation, instance)) {
                    offspring.push_back(child);
                }
            }
            
            // Create second offspring (reverse roles)
            std::vector<int> offspring2Genes(genes1.size());
            std::vector<bool> used2(genes1.size(), false);
            
            for (int i = start; i <= end; i++) {
                offspring2Genes[i] = genes2[i];
                used2[genes2[i]] = true;
            }
            
            j = 0;
            bool validOffspring2 = true;
            for (int i = 0; i < static_cast<int>(genes1.size()); i++) {
                if (i >= start && i <= end) continue;
                
                // Find next unused gene from parent1
                while (j < static_cast<int>(genes1.size()) && used2[genes1[j]]) j++;
                
                if (j < static_cast<int>(genes1.size())) {
                    offspring2Genes[i] = genes1[j++];
                } else {
                    // If we run out of genes, find any unused valid index
                    bool foundUnused = false;
                    for (size_t val = 0; val < instance.getSpectrum().size(); val++) {
                        if (!used2[val]) {
                            offspring2Genes[i] = val;
                            used2[val] = true;
                            foundUnused = true;
                            break;
                        }
                    }
                    if (!foundUnused) {
                        validOffspring2 = false;
                        break;
                    }
                }
            }
            
            if (validOffspring2) {
                if (auto child = createOffspring(offspring2Genes, representation, instance)) {
                    offspring.push_back(child);
                }
            }
            
            maxAttempts--;
        }
        
        // If we couldn't create any valid offspring, return copies of parents
        if (offspring.empty()) {
            LOG_WARNING("Failed to create valid offspring, returning parents");
            return parents;
        }
        
        return offspring;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in order crossover: " + std::string(e.what()));
        return parents;  // Return parents instead of empty vector
    }
}

// ========================================
// =       EDGE RECOMBINATION (ERX)      =
// ========================================
std::vector<std::shared_ptr<Individual>> 
EdgeRecombination::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
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
                    if (!used[neighbor]) {
                        int unusedNeighborCount = 0;
                        for (int n : edgeTable[neighbor]) {
                            if (!used[n]) unusedNeighborCount++;
                        }
                        if (unusedNeighborCount < minNeighbors) {
                            minNeighbors = unusedNeighborCount;
                            next = neighbor;
                        }
                    }
                }
                
                // If no unused neighbors found, pick random unused gene
                if (next == -1) {
                    std::vector<int> unusedGenes;
                    for (size_t i = 0; i < genes1.size(); i++) {
                        if (!used[i]) unusedGenes.push_back(i);
                    }
                    if (!unusedGenes.empty()) {
                        next = unusedGenes[rng.getRandomInt(0, static_cast<int>(unusedGenes.size() - 1))];
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

// ========================================
// =       PMX CROSSOVER                 =
// ========================================
std::vector<std::shared_ptr<Individual>> 
PMXCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
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
        
        // Create mapping from parent2 to parent1 values
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
        
        // Create second offspring
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