//
// Created by konrad_guest on 07/01/2025.
//
#include "metaheuristics/crossover.h"
#include <random>

// ============ OnePointCrossover ============
std::vector<void*> 
OnePointCrossover::crossover(const std::vector<void*>        &parents,
                             const DNAInstance               &instance,
                             std::shared_ptr<IRepresentation> representation)
{
    // Example stub: no actual crossover, just copy parents
    std::vector<void*> offspring = parents;
    return offspring;
}

// ----------------------------------------------------------------------
//                           Order Crossover (OX)
// ----------------------------------------------------------------------
std::vector<void*> 
OrderCrossover::crossover(const std::vector<void*> &parents,
                          const DNAInstance &instance,
                          std::shared_ptr<IRepresentation> representation)
{
    // Zakładamy, że parents to [p0, p1, p2, p3, p4, p5, ...] w parach
    // Będziemy generować tyle dzieci, by zastąpić populację (lub tyle samo).
    std::vector<void*> offspring;
    offspring.reserve(parents.size()); // jeżeli 1 potomek na parę

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto *parent1 = static_cast<std::vector<int>*>(parents[i]);
        auto *parent2 = static_cast<std::vector<int>*>(parents[i+1]);
        if(!parent1 || !parent2) continue;
        if(parent1->size() != parent2->size()) {
            std::cerr << "OrderCrossover: mismatch sizes" << std::endl;
            continue;
        }

        int size = (int)parent1->size();

        // (Tworzymy 1 dziecko – OX)
        std::vector<int> child(size, -1);
        std::vector<bool> used(size, false);
        int start = rand() % size;
        int end   = rand() % size;
        if(start > end) std::swap(start, end);

        // Kopiujemy [start..end] z parent1
        for(int j=start; j<=end; j++){
            child[j] = (*parent1)[j];
            used[(*parent1)[j]] = true;
        }

        // Uzupełniamy resztę z parent2
        int curr = (end+1) % size;
        int p2pos= (end+1) % size;
        while(curr != start){
            while( used[(*parent2)[p2pos]] ) {
                p2pos = (p2pos+1) % size;
            }
            child[curr] = (*parent2)[p2pos];
            used[(*parent2)[p2pos]] = true;
            curr = (curr+1) % size;
            p2pos= (p2pos+1) % size;
        }

        // dodajemy do offspring
        auto *childPtr = new std::vector<int>(child);
        offspring.push_back(childPtr);

        // [Opcjonalnie: stwórzmy 2. dziecko z parent2->parent1, by zachować populację w 1:1]
        // ...
    }

    return offspring;
}

// ----------------------------------------------------------------------
//                  Edge Recombination Crossover (ERX)
// ----------------------------------------------------------------------
std::vector<void*> 
EdgeRecombination::crossover(const std::vector<void*>        &parents,
                             const DNAInstance               &instance,
                             std::shared_ptr<IRepresentation> representation)
{
    if (parents.size() < 2) return parents;
    
    auto* parent1 = static_cast<std::vector<int>*>(parents[0]);
    auto* parent2 = static_cast<std::vector<int>*>(parents[1]);
    
    int size = static_cast<int>(parent1->size());
    std::vector<std::vector<int>> edgeMap(size);
    std::vector<int> child;
    child.reserve(size);
    std::vector<bool> used(size, false);
    
    // Helper lambda to add edges from a parent
    auto addEdges = [&](const std::vector<int> &p) {
        for (int i = 0; i < size; i++) {
            int curr = p[i];
            int next = p[(i + 1) % size];
            int prev = p[(i - 1 + size) % size];
            // add edges from curr -> prev and curr -> next
            edgeMap[curr].push_back(prev);
            edgeMap[curr].push_back(next);
        }
    };
    
    // Build edgeMap from both parents
    addEdges(*parent1);
    addEdges(*parent2);
    
    // Remove duplicates in each adjacency list
    for (auto &edges : edgeMap) {
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    }
    
    // Start with parent1[0] 
    int current = (*parent1)[0];
    child.push_back(current);
    used[current] = true;
    
    while ((int)child.size() < size) {
        // gather unused neighbors of current
        std::vector<int> neighbors;
        for (int nb : edgeMap[current]) {
            if (!used[nb]) {
                neighbors.push_back(nb);
            }
        }
        
        if (neighbors.empty()) {
            // pick random unused gene
            std::vector<int> unusedGenes;
            for (int i = 0; i < size; i++) {
                if (!used[i]) {
                    unusedGenes.push_back(i);
                }
            }
            if (!unusedGenes.empty()) {
                current = unusedGenes[ rand() % unusedGenes.size() ];
            } else {
                // Should never happen if we haven't used all genes
                break;
            }
        } else {
            // choose a neighbor randomly
            current = neighbors[ rand() % neighbors.size() ];
        }
        
        child.push_back(current);
        used[current] = true;
    }
    
    std::vector<void*> offspring;
    offspring.push_back( new std::vector<int>(child) );
    return offspring;
}

// ----------------------------------------------------------------------
//                 Partially Mapped Crossover (PMX)
// ----------------------------------------------------------------------
std::vector<void*> 
PMXCrossover::crossover(const std::vector<void*>        &parents,
                        const DNAInstance               &instance,
                        std::shared_ptr<IRepresentation> representation)
{
    if (parents.size() < 2) return parents;
    
    auto* parent1 = static_cast<std::vector<int>*>(parents[0]);
    auto* parent2 = static_cast<std::vector<int>*>(parents[1]);
    
    int size = static_cast<int>(parent1->size());
    std::vector<int> child(*parent1);  // copy parent1
    std::vector<int> position(size);
    
    // Build position table for parent2
    // position[val] = index of val in parent2
    for (int i = 0; i < size; i++) {
        position[ (*parent2)[i] ] = i;
    }
    
    // pick random crossover range
    int start = rand() % size;
    int end   = rand() % size;
    if (start > end) std::swap(start, end);
    
    // PMX: swap conflicting genes in the slice [start..end]
    for (int i = start; i <= end; i++) {
        int val1 = child[i];
        int val2 = (*parent2)[i];
        
        if (val1 != val2) {
            // positions of val1, val2
            int pos1 = i;
            int pos2 = position[val1];
            
            // swap them in child
            std::swap(child[pos1], child[pos2]);
            
            // update position
            position[val1] = pos2;
            position[val2] = pos1;
        }
    }
    
    std::vector<void*> offspring;
    offspring.push_back( new std::vector<int>(child) );
    return offspring;
}



// // ============ TwoPointCrossover ============
// std::vector<std::vector<double>>
// TwoPointCrossover::crossover(const std::vector<std::vector<double>> &parents)
// {
//     // ...
//     return {};
// }
//
// // ============ UniformCrossover ============
// std::vector<std::vector<double>>
// UniformCrossover::crossover(const std::vector<std::vector<double>> &parents)
// {
//     // ...
//     return {};
// }
//
// // ============ ArithmeticCrossover ============
// std::vector<std::vector<double>>
// ArithmeticCrossover::crossover(const std::vector<std::vector<double>> &parents)
// {
//     // ...
//     return {};
// }