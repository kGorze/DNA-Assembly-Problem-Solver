#include "metaheuristics/crossover_impl.h"
#include "utils/logger.h"
#include <algorithm>
#include <random>
#include <unordered_set>
#include <limits>

std::vector<std::shared_ptr<Individual>> OnePointCrossover::crossover(
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
    
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < parents.size() - 1; i += 2) {
        if (dis(gen) < m_crossoverRate) {
            auto [child1, child2] = performCrossover(parents[i], parents[i + 1], instance, representation);
            if (child1) offspring.push_back(child1);
            if (child2) offspring.push_back(child2);
        } else {
            offspring.push_back(parents[i]);
            offspring.push_back(parents[i + 1]);
        }
    }
    
    // Handle odd number of parents
    if (parents.size() % 2 == 1) {
        offspring.push_back(parents.back());
    }
    
    return offspring;
}

std::pair<std::shared_ptr<Individual>, std::shared_ptr<Individual>>
OnePointCrossover::performCrossover(
    const std::shared_ptr<Individual>& parent1,
    const std::shared_ptr<Individual>& parent2,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (!parent1 || !parent2) {
        LOG_ERROR("Null parent in crossover");
        return {nullptr, nullptr};
    }
    
    const auto& genes1 = parent1->getGenes();
    const auto& genes2 = parent2->getGenes();
    
    if (genes1.empty() || genes2.empty()) {
        LOG_ERROR("Empty parent in crossover");
        return {nullptr, nullptr};
    }
    
    std::uniform_int_distribution<> dis(1, std::min(genes1.size(), genes2.size()) - 1);
    size_t crossPoint = dis(gen);
    
    std::vector<int> offspring1;
    std::vector<int> offspring2;
    offspring1.reserve(genes1.size());
    offspring2.reserve(genes2.size());
    
    // First child gets first part from parent1, second part from parent2
    offspring1.insert(offspring1.end(), genes1.begin(), genes1.begin() + crossPoint);
    offspring1.insert(offspring1.end(), genes2.begin() + crossPoint, genes2.end());
    
    // Second child gets first part from parent2, second part from parent1
    offspring2.insert(offspring2.end(), genes2.begin(), genes2.begin() + crossPoint);
    offspring2.insert(offspring2.end(), genes1.begin() + crossPoint, genes1.end());
    
    auto child1 = std::make_shared<Individual>(offspring1);
    auto child2 = std::make_shared<Individual>(offspring2);
    
    if (!representation->isValid(child1) || !representation->isValid(child2)) {
        LOG_ERROR("Invalid offspring produced in crossover");
        return {nullptr, nullptr};
    }
    
    return {child1, child2};
}

std::vector<std::shared_ptr<Individual>> OrderCrossover::crossover(
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
        
        auto child1Genes = performOrderCrossover(genes1, genes2, size);
        auto child2Genes = performOrderCrossover(genes2, genes1, size);
        
        auto child1 = std::make_shared<Individual>(child1Genes);
        auto child2 = std::make_shared<Individual>(child2Genes);
        
        if (representation->isValid(child1)) {
            offspring.push_back(child1);
        }
        if (representation->isValid(child2)) {
            offspring.push_back(child2);
        }
    }
    
    // Handle odd number of parents
    if (parents.size() % 2 == 1) {
        offspring.push_back(parents.back());
    }
    
    return offspring;
}

std::vector<int> OrderCrossover::performOrderCrossover(
    const std::vector<int>& parent1,
    const std::vector<int>& parent2,
    size_t size)
{
    std::vector<int> child(size, -1);
    std::vector<bool> used(size, false);
    
    // Select random subsequence from parent1
    std::uniform_int_distribution<> dis(0, size - 1);
    size_t start = dis(gen);
    size_t length = dis(gen) + 1;
    
    // Copy subsequence from parent1
    for (size_t i = 0; i < length; ++i) {
        size_t pos = (start + i) % size;
        child[pos] = parent1[pos];
        used[parent1[pos]] = true;
    }
    
    // Fill remaining positions with genes from parent2 in order
    size_t childPos = (start + length) % size;
    for (size_t i = 0; i < size; ++i) {
        size_t parentPos = (start + length + i) % size;
        if (!used[parent2[parentPos]]) {
            child[childPos] = parent2[parentPos];
            childPos = (childPos + 1) % size;
        }
    }
    
    return child;
}

EdgeRecombination::EdgeTable::EdgeTable(
    const std::vector<std::shared_ptr<Individual>>& parents)
{
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
        
        // Create edge table from both parents
        EdgeTable edgeTable({parent1, parent2});
        
        // Create offspring
        std::vector<int> childGenes;
        childGenes.reserve(genes1.size());
        
        // Start with a random node from parent1
        std::uniform_int_distribution<> dis(0, genes1.size() - 1);
        int currentNode = genes1[dis(gen)];
        childGenes.push_back(currentNode);
        edgeTable.removeNode(currentNode);
        
        // Build the rest of the path
        while (childGenes.size() < genes1.size()) {
            auto neighbors = edgeTable.getNeighbors(currentNode);
            
            if (neighbors.empty()) {
                // If no neighbors, choose random remaining node
                std::vector<int> remainingNodes;
                for (const auto& node : genes1) {
                    if (edgeTable.hasNode(node)) {
                        remainingNodes.push_back(node);
                    }
                }
                
                if (remainingNodes.empty()) break;
                
                std::uniform_int_distribution<> dis(0, remainingNodes.size() - 1);
                currentNode = remainingNodes[dis(gen)];
            } else {
                // Choose neighbor with fewest remaining neighbors
                int bestNode = neighbors[0];
                size_t minNeighbors = std::numeric_limits<size_t>::max();
                
                for (int node : neighbors) {
                    size_t numNeighbors = edgeTable.getNeighbors(node).size();
                    if (numNeighbors < minNeighbors) {
                        minNeighbors = numNeighbors;
                        bestNode = node;
                    }
                }
                
                currentNode = bestNode;
            }
            
            childGenes.push_back(currentNode);
            edgeTable.removeNode(currentNode);
        }
        
        auto child = std::make_shared<Individual>(childGenes);
        
        if (representation->isValid(child)) {
            offspring.push_back(child);
        }
    }
    
    // Handle odd number of parents
    if (parents.size() % 2 == 1) {
        offspring.push_back(parents.back());
    }
    
    return offspring;
} 