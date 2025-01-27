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
    
    // Basic validation
    if (parents.size() < 2) {
        LOG_WARNING("OrderCrossover: Not enough parents for crossover");
        return parents;
    }
    
    if (!parents[0] || !parents[1]) {
        LOG_WARNING("OrderCrossover: Null parent received");
        return parents;
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    
    if (parent1Genes.empty() || parent2Genes.empty()) {
        LOG_WARNING("OrderCrossover: Empty genes in parent");
        return parents;
    }
    
    // Create offspring
    std::vector<std::shared_ptr<Individual>> offspring;
    
    // Get size from parent genes
    size_t size = parent1Genes.size();
    
    // Try to create valid offspring multiple times
    const int maxAttempts = 5;
    int attempts = 0;
    bool success = false;
    
    while (attempts < maxAttempts && !success) {
        try {
            // Perform crossover operation with size parameter
            std::vector<int> child1Genes = performOrderCrossover(parent1Genes, parent2Genes, size);
            std::vector<int> child2Genes = performOrderCrossover(parent2Genes, parent1Genes, size);
            
            // Basic validation of gene indices
            bool genesValid = true;
            for (int gene : child1Genes) {
                if (gene < 0 || gene >= static_cast<int>(instance.getSpectrum().size())) {
                    genesValid = false;
                    break;
                }
            }
            
            if (!genesValid) {
                LOG_WARNING("OrderCrossover: Invalid gene indices in offspring");
                attempts++;
                continue;
            }
            
            // Create offspring individuals
            auto child1 = std::make_shared<Individual>(child1Genes);
            auto child2 = std::make_shared<Individual>(child2Genes);
            
            // Validate offspring
            if (representation->isValid(child1, instance) && representation->isValid(child2, instance)) {
                offspring.push_back(child1);
                offspring.push_back(child2);
                success = true;
            } else {
                LOG_WARNING("OrderCrossover: Invalid offspring produced");
                attempts++;
            }
            
        } catch (const std::exception& e) {
            LOG_WARNING("OrderCrossover: Exception during crossover: " + std::string(e.what()));
            attempts++;
        }
    }
    
    // If no valid offspring were produced, return copies of parents
    if (offspring.empty()) {
        LOG_WARNING("OrderCrossover: Failed to produce valid offspring after " + std::to_string(maxAttempts) + " attempts");
        return parents;
    }
    
    return offspring;
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
        int stuckCounter = 0;
        while (childGenes.size() < parent1.size() && stuckCounter < 100) {
            const auto neighbors = edgeTable.getNeighbors(current);
            bool foundNext = false;
            
            // First try: Find unused neighbor with fewest remaining neighbors
            if (!neighbors.empty()) {
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
                        foundNext = true;
                    }
                }
                
                if (foundNext && next != -1) {
                    childGenes.push_back(next);
                    used.insert(next);
                    current = next;
                    edgeTable.removeNode(current);
                    stuckCounter = 0;
                    continue;
                }
            }
            
            // Second try: Use any unused gene from parent1 or parent2
            std::vector<int> unusedGenes;
            for (int gene : parent1) {
                if (used.find(gene) == used.end()) {
                    unusedGenes.push_back(gene);
                }
            }
            for (int gene : parent2) {
                if (used.find(gene) == used.end() && 
                    std::find(unusedGenes.begin(), unusedGenes.end(), gene) == unusedGenes.end()) {
                    unusedGenes.push_back(gene);
                }
            }
            
            if (!unusedGenes.empty()) {
                std::uniform_int_distribution<size_t> unusedDist(0, unusedGenes.size() - 1);
                int next = unusedGenes[unusedDist(gen)];
                childGenes.push_back(next);
                used.insert(next);
                current = next;
                edgeTable.removeNode(current);
                stuckCounter = 0;
            } else {
                // Last resort: Use any unused value from 0 to size-1
                for (size_t i = 0; i < parent1.size(); i++) {
                    if (used.find(i) == used.end()) {
                        childGenes.push_back(i);
                        used.insert(i);
                        current = i;
                        edgeTable.removeNode(current);
                        stuckCounter = 0;
                        break;
                    }
                }
            }
            
            stuckCounter++;
        }
        
        // If we got stuck, fill remaining positions with unused values
        if (childGenes.size() < parent1.size()) {
            for (size_t i = 0; i < parent1.size() && childGenes.size() < parent1.size(); i++) {
                if (used.find(i) == used.end()) {
                    childGenes.push_back(i);
                    used.insert(i);
                }
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
    
    // Continue until we reach point1 again
    while (curr != point1) {
        bool foundUnused = false;
        
        // Try one complete cycle through parent2 starting from p2pos
        for (size_t i = 0; i < size; ++i) {
            size_t checkPos = (p2pos + i) % size;
            if (used.find(parent2[checkPos]) == used.end()) {
                child[curr] = parent2[checkPos];
                used.insert(parent2[checkPos]);
                p2pos = (checkPos + 1) % size;
                foundUnused = true;
                break;
            }
        }
        
        // If no unused gene found in parent2, find any unused value
        if (!foundUnused) {
            for (size_t val = 0; val < size; ++val) {
                if (used.find(val) == used.end()) {
                    child[curr] = val;
                    used.insert(val);
                    break;
                }
            }
        }
        
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
        
        // Track used values for each child
        std::unordered_set<int> used1, used2;
        
        // Copy the mapping section
        for (size_t j = point1; j <= point2; ++j) {
            child1[j] = genes2[j];
            child2[j] = genes1[j];
            used1.insert(child1[j]);
            used2.insert(child2[j]);
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
            int safetyCounter1 = 0;
            std::unordered_set<int> visited1;
            
            while (mapping1.find(value1) != mapping1.end() && safetyCounter1 < size) {
                if (visited1.find(value1) != visited1.end()) {
                    // We've found a cycle, break it by using any unused value
                    for (size_t val = 0; val < size; ++val) {
                        if (used1.find(val) == used1.end()) {
                            value1 = val;
                            break;
                        }
                    }
                    break;
                }
                visited1.insert(value1);
                value1 = mapping1[value1];
                safetyCounter1++;
            }
            
            // If we still haven't found a valid value, use any unused value
            if (used1.find(value1) != used1.end()) {
                for (size_t val = 0; val < size; ++val) {
                    if (used1.find(val) == used1.end()) {
                        value1 = val;
                        break;
                    }
                }
            }
            
            child1[j] = value1;
            used1.insert(value1);
            
            // For child2
            int value2 = genes2[j];
            int safetyCounter2 = 0;
            std::unordered_set<int> visited2;
            
            while (mapping2.find(value2) != mapping2.end() && safetyCounter2 < size) {
                if (visited2.find(value2) != visited2.end()) {
                    // We've found a cycle, break it by using any unused value
                    for (size_t val = 0; val < size; ++val) {
                        if (used2.find(val) == used2.end()) {
                            value2 = val;
                            break;
                        }
                    }
                    break;
                }
                visited2.insert(value2);
                value2 = mapping2[value2];
                safetyCounter2++;
            }
            
            // If we still haven't found a valid value, use any unused value
            if (used2.find(value2) != used2.end()) {
                for (size_t val = 0; val < size; ++val) {
                    if (used2.find(val) == used2.end()) {
                        value2 = val;
                        break;
                    }
                }
            }
            
            child2[j] = value2;
            used2.insert(value2);
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