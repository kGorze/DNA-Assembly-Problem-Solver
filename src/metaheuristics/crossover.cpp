//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/crossover.h"
#include <random>
#include <iostream>
#include <algorithm>

// ============ OnePointCrossover ============

std::vector<std::shared_ptr<std::vector<int>>> 
OnePointCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    // Przykładowa, prosta implementacja – kopiuje rodziców
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());
    for (auto &p : parents) {
        if (p) {
            auto child = std::make_shared<std::vector<int>>(*p);
            offspring.push_back(child);
        }
    }
    return offspring;
}

// ----------------------------------------------------------------------
// Order Crossover (OX)
// ----------------------------------------------------------------------
std::vector<std::shared_ptr<std::vector<int>>> 
OrderCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());

    std::random_device rd;
    std::mt19937 gen(rd());

    // Krzyżujemy parami
    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto parent1 = parents[i];
        auto parent2 = parents[i+1];
        // DODANE ZABEZPIECZENIA
        if (!parent1 || !parent2 || parent1->empty() || parent2->empty()) {
            continue;
        }
        if (parent1->size() != parent2->size()) {
            std::cerr << "OrderCrossover: mismatch in parent sizes" << std::endl;
            continue;
        }

        int size = static_cast<int>(parent1->size());
        std::uniform_int_distribution<int> distIndex(0, size - 1);
        int start = distIndex(gen);
        int end   = distIndex(gen);
        if (start > end) std::swap(start, end);

        // Tworzymy dziecko
        std::vector<int> child(size, -1);
        std::vector<bool> used(size, false);

        for (int j = start; j <= end; j++) {
            child[j] = (*parent1)[j];
            used[(*parent1)[j]] = true;
        }

        int curr = (end + 1) % size;
        int p2pos = (end + 1) % size;
        while (curr != start) {
            while (used[(*parent2)[p2pos]]) {
                p2pos = (p2pos + 1) % size;
            }
            child[curr] = (*parent2)[p2pos];
            used[(*parent2)[p2pos]] = true;
            curr = (curr + 1) % size;
            p2pos = (p2pos + 1) % size;
        }

        offspring.push_back(std::make_shared<std::vector<int>>(child));
    }

    return offspring;
}

// ----------------------------------------------------------------------
// Edge Recombination Crossover (ERX)
// ----------------------------------------------------------------------
std::vector<std::shared_ptr<std::vector<int>>> 
EdgeRecombination::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    // DODANE ZABEZPIECZENIA
    if (parents.size() < 2) {
        return parents;
    }

    auto parent1 = parents[0];
    auto parent2 = parents[1];
    if (!parent1 || !parent2 || parent1->empty() || parent2->empty()) {
        return parents;
    }
    if (parent1->size() != parent2->size()) {
        std::cerr << "EdgeRecombination: mismatch parent sizes." << std::endl;
        return {};
    }

    int size = static_cast<int>(parent1->size());
    std::vector<std::vector<int>> edgeMap(size);
    std::vector<int> child;
    child.reserve(size);
    std::vector<bool> used(size, false);

    auto addEdges = [&](const std::shared_ptr<std::vector<int>>& p) {
        for (int i = 0; i < size; i++) {
            int curr = (*p)[i];
            int next = (*p)[(i + 1) % size];
            int prev = (*p)[(i - 1 + size) % size];
            edgeMap[curr].push_back(prev);
            edgeMap[curr].push_back(next);
        }
    };

    addEdges(parent1);
    addEdges(parent2);

    for (auto &edges : edgeMap) {
        std::sort(edges.begin(), edges.end());
        edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    }

    std::random_device rd;
    std::mt19937 gen(rd());

    int current = (*parent1)[0];
    child.push_back(current);
    used[current] = true;

    while ((int)child.size() < size) {
        std::vector<int> neighbors;
        for (int nb : edgeMap[current]) {
            if (!used[nb]) {
                neighbors.push_back(nb);
            }
        }

        if (neighbors.empty()) {
            std::vector<int> unusedGenes;
            unusedGenes.reserve(size);
            for (int i = 0; i < size; i++) {
                if (!used[i]) {
                    unusedGenes.push_back(i);
                }
            }
            if (!unusedGenes.empty()) {
                std::uniform_int_distribution<int> distU(0, (int)unusedGenes.size() - 1);
                current = unusedGenes[distU(gen)];
            } else {
                break;
            }
        } else {
            std::uniform_int_distribution<int> distN(0, (int)neighbors.size() - 1);
            current = neighbors[distN(gen)];
        }

        child.push_back(current);
        used[current] = true;
    }

    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.push_back(std::make_shared<std::vector<int>>(child));
    return offspring;
}

// ----------------------------------------------------------------------
// Partially Mapped Crossover (PMX)
// ----------------------------------------------------------------------
std::vector<std::shared_ptr<std::vector<int>>> 
PMXCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    // DODANE ZABEZPIECZENIA
    if (parents.size() < 2) {
        return {};
    }

    auto parent1 = parents[0];
    auto parent2 = parents[1];
    if (!parent1 || !parent2 || parent1->empty() || parent2->empty()) {
        return {};
    }
    if (parent1->size() != parent2->size()) {
        return {};
    }

    int size = static_cast<int>(parent1->size());
    std::vector<int> child(*parent1);
    std::vector<int> position(size);

    for (int i = 0; i < size; i++) {
        position[(*parent2)[i]] = i;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distIndex(0, size - 1);

    int start = distIndex(gen);
    int end   = distIndex(gen);
    if (start > end) std::swap(start, end);

    for (int i = start; i <= end; i++) {
        int val1 = child[i];
        int val2 = (*parent2)[i];
        if (val1 != val2) {
            int pos1 = i;
            int pos2 = position[val1];
            std::swap(child[pos1], child[pos2]);
            position[val1] = pos2;
            position[val2] = pos1;
        }
    }

    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.push_back(std::make_shared<std::vector<int>>(child));
    return offspring;
}
