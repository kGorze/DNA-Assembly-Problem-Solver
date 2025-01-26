#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/individual.h"
#include <algorithm>
#include <random>
#include <unordered_set>
#include <limits>
#include <unordered_map>

// Forward declaration of helper function
std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size);

std::vector<std::shared_ptr<Individual>> OrderCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !representation) {
        LOG_WARNING("Invalid input for crossover");
        return parents;
    }

    const auto& parent1 = parents[0];
    const auto& parent2 = parents[1];
    
    if (!parent1 || !parent2) {
        LOG_WARNING("Null parents in crossover");
        return parents;
    }

    const auto& p1Genes = parent1->getGenes();
    const auto& p2Genes = parent2->getGenes();
    
    if (p1Genes.empty() || p2Genes.empty()) {
        LOG_WARNING("Empty genes in crossover");
        return parents;
    }

    try {
        // Calculate allowed mismatches based on instance parameters
        const int deltaK = instance.getDeltaK();
        const int totalErrors = instance.getLNeg() + instance.getLPoz();
        const int extraMismatches = (deltaK == 0 && totalErrors > 0) ? 
            std::min(2, totalErrors / 10) : 0;
        const int allowedMismatches = deltaK + extraMismatches;

        std::random_device rd;
        std::mt19937 gen(rd());
        
        // Select crossover points
        std::uniform_int_distribution<size_t> dis(0, p1Genes.size() - 1);
        size_t point1 = dis(gen);
        size_t point2 = dis(gen);
        
        if (point1 > point2) {
            std::swap(point1, point2);
        }
        
        // Create offspring
        std::vector<int> offspring1Genes(p1Genes.size(), -1);
        std::vector<int> offspring2Genes(p2Genes.size(), -1);
        
        // Copy selected segments
        std::copy(p1Genes.begin() + point1, p1Genes.begin() + point2 + 1, 
                 offspring1Genes.begin() + point1);
        std::copy(p2Genes.begin() + point1, p2Genes.begin() + point2 + 1, 
                 offspring2Genes.begin() + point1);
        
        // Fill remaining positions with relaxed validation
        auto fillOffspring = [&](const std::vector<int>& sourceGenes, 
                               std::vector<int>& targetGenes) {
            std::vector<bool> used(instance.getSpectrum().size(), false);
            
            // Mark used genes in the crossover segment
            for (size_t i = point1; i <= point2; i++) {
                if (targetGenes[i] >= 0) {
                    used[targetGenes[i]] = true;
                }
            }
            
            // Fill positions before crossover point
            size_t j = 0;
            for (size_t i = 0; i < point1; i++) {
                while (j < sourceGenes.size() && used[sourceGenes[j]]) j++;
                if (j < sourceGenes.size()) {
                    targetGenes[i] = sourceGenes[j++];
                    used[targetGenes[i]] = true;
                }
            }
            
            // Fill positions after crossover point
            for (size_t i = point2 + 1; i < targetGenes.size(); i++) {
                while (j < sourceGenes.size() && used[sourceGenes[j]]) j++;
                if (j < sourceGenes.size()) {
                    targetGenes[i] = sourceGenes[j++];
                    used[targetGenes[i]] = true;
                }
            }
        };
        
        // Create both offspring
        fillOffspring(p2Genes, offspring1Genes);
        fillOffspring(p1Genes, offspring2Genes);
        
        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);
        
        // Create and add offspring without strict validation
        auto child1 = std::make_shared<Individual>(offspring1Genes);
        auto child2 = std::make_shared<Individual>(offspring2Genes);
        
        // Always add offspring, let fitness function handle penalties
        offspring.push_back(child1);
        offspring.push_back(child2);
        
        return offspring;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in crossover: " + std::string(e.what()));
        return parents;
    }
}

// EdgeTable implementation
EdgeRecombination::EdgeTable::EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents) {
    if (parents.empty()) return;
    
    // Initialize edge lists for each node
    for (const auto& parent : parents) {
        if (!parent) continue;
        
        const auto& genes = parent->getGenes();
        if (genes.empty()) continue;
        
        for (size_t i = 0; i < genes.size(); ++i) {
            int current = genes[i];
            int prev = i > 0 ? genes[i - 1] : genes.back();
            int next = i < genes.size() - 1 ? genes[i + 1] : genes.front();
            
            edges[current].push_back(prev);
            edges[current].push_back(next);
        }
    }
    
    // Remove duplicates from edge lists
    for (auto& [node, neighbors] : edges) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
}

std::vector<int> EdgeRecombination::EdgeTable::getNeighbors(int node) const {
    auto it = edges.find(node);
    return it != edges.end() ? it->second : std::vector<int>();
}

void EdgeRecombination::EdgeTable::removeNode(int node) {
    // Remove node from all neighbor lists
    for (auto& [_, neighbors] : edges) {
        auto it = std::remove(neighbors.begin(), neighbors.end(), node);
        neighbors.erase(it, neighbors.end());
    }
    // Remove node's edge list
    edges.erase(node);
}

bool EdgeRecombination::EdgeTable::hasNode(int node) const {
    return edges.find(node) != edges.end();
}

std::vector<std::shared_ptr<Individual>> EdgeRecombination::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    if (parents.size() < 2) {
        LOG_WARNING("Insufficient parents for crossover");
        return parents;
    }

    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);

    const auto& parent1 = parents[0]->getGenes();
    const auto& parent2 = parents[1]->getGenes();
    
    if (parent1.empty() || parent2.empty()) {
        LOG_WARNING("Empty parent genes in crossover");
        return parents;
    }

    // Build edge table
    EdgeTable edgeTable(parents);
    
    // Create two offspring
    for (int k = 0; k < 2; k++) {
        std::vector<int> childGenes;
        childGenes.reserve(parent1.size());
        std::unordered_set<int> used;
        
        // Start with a random gene from first parent
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(0, parent1.size() - 1);
        int current = parent1[dist(gen)];
        childGenes.push_back(current);
        used.insert(current);
        
        // Build the rest of the path
        while (childGenes.size() < parent1.size()) {
            const auto neighbors = edgeTable.getNeighbors(current);
            
            // Find unused neighbor with fewest remaining neighbors
            int next = -1;
            int minNeighbors = std::numeric_limits<int>::max();
            
            for (int neighbor : neighbors) {
                if (used.find(neighbor) != used.end()) continue;
                
                auto neighborEdges = edgeTable.getNeighbors(neighbor);
                int remainingNeighbors = 0;
                for (int n : neighborEdges) {
                    if (used.find(n) == used.end()) remainingNeighbors++;
                }
                
                if (remainingNeighbors < minNeighbors) {
                    minNeighbors = remainingNeighbors;
                    next = neighbor;
                }
            }
            
            // If no unused neighbors, select random unused gene
            if (next == -1) {
                std::vector<int> unusedGenes;
                for (const auto& gene : parent1) {
                    if (used.find(gene) == used.end()) {
                        unusedGenes.push_back(gene);
                    }
                }
                
                if (!unusedGenes.empty()) {
                    std::uniform_int_distribution<size_t> unusedDist(0, unusedGenes.size() - 1);
                    next = unusedGenes[unusedDist(gen)];
                }
            }
            
            if (next != -1) {
                childGenes.push_back(next);
                used.insert(next);
                current = next;
                edgeTable.removeNode(current);
            }
        }
        
        auto child = std::make_shared<Individual>(childGenes);
        if (representation->isValid(child, instance)) {
            offspring.push_back(child);
        }
    }
    
    // If we couldn't create valid offspring, return copies of parents
    if (offspring.empty()) {
        offspring.push_back(parents[0]);
        offspring.push_back(parents[1]);
    }
    
    return offspring;
}

std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, size - 1);
    
    // Select two random crossover points
    size_t point1 = dist(gen);
    size_t point2 = dist(gen);
    if (point1 > point2) std::swap(point1, point2);
    
    // Create child with same size as parents
    std::vector<int> child(size, -1);
    std::unordered_set<int> used;
    
    // Copy segment between crossover points from parent1
    for (size_t i = point1; i <= point2; ++i) {
        child[i] = parent1[i];
        used.insert(parent1[i]);
    }
    
    // Fill remaining positions with genes from parent2 in order
    size_t curr = (point2 + 1) % size;
    size_t p2pos = (point2 + 1) % size;
    
    while (curr != point1) {
        // Find next unused gene from parent2
        while (used.find(parent2[p2pos]) != used.end()) {
            p2pos = (p2pos + 1) % size;
        }
        
        child[curr] = parent2[p2pos];
        used.insert(parent2[p2pos]);
        curr = (curr + 1) % size;
    }
    
    return child;
}

std::vector<std::shared_ptr<Individual>> PMXCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (parents.size() < 2) {
        LOG_ERROR("Not enough parents for crossover");
        return {};
    }
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(parents.size());
    
    for (size_t i = 0; i < parents.size() - 1; i += 2) {
        const auto& parent1 = parents[i];
        const auto& parent2 = parents[i + 1];
        
        if (!parent1 || !parent2) {
            LOG_ERROR("Null parent in crossover");
            continue;
        }
        
        const auto& genes1 = parent1->getGenes();
        const auto& genes2 = parent2->getGenes();
        
        if (genes1.empty() || genes2.empty()) {
            LOG_ERROR("Empty parent in crossover");
            continue;
        }
        
        size_t size = std::min(genes1.size(), genes2.size());
        
        // Generate two random crossover points
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(0, size - 1);
        size_t point1 = dist(gen);
        size_t point2 = dist(gen);
        if (point1 > point2) std::swap(point1, point2);
        
        // Create offspring
        std::vector<int> child1(size, -1);
        std::vector<int> child2(size, -1);
        
        // Copy the mapping section
        for (size_t j = point1; j <= point2; ++j) {
            child1[j] = genes2[j];
            child2[j] = genes1[j];
        }
        
        // Create mapping between values in the mapping section
        std::unordered_map<int, int> mapping1, mapping2;
        for (size_t j = point1; j <= point2; ++j) {
            mapping1[genes1[j]] = genes2[j];
            mapping2[genes2[j]] = genes1[j];
        }
        
        // Fill in remaining positions
        for (size_t j = 0; j < size; ++j) {
            if (j >= point1 && j <= point2) continue;
            
            // For child1
            int value1 = genes1[j];
            while (mapping1.find(value1) != mapping1.end()) {
                value1 = mapping1[value1];
            }
            child1[j] = value1;
            
            // For child2
            int value2 = genes2[j];
            while (mapping2.find(value2) != mapping2.end()) {
                value2 = mapping2[value2];
            }
            child2[j] = value2;
        }
        
        auto offspring1 = std::make_shared<Individual>(child1);
        auto offspring2 = std::make_shared<Individual>(child2);
        
        if (representation->isValid(offspring1, instance)) {
            offspring.push_back(offspring1);
        }
        if (representation->isValid(offspring2, instance)) {
            offspring.push_back(offspring2);
        }
    }
    
    // Handle odd number of parents
    if (parents.size() % 2 == 1) {
        offspring.push_back(parents.back());
    }
    
    return offspring;
} 