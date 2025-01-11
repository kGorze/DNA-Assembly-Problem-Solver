//
// Created by konrad_guest on 07/01/2025.
// SMART

#include "metaheuristics/crossover.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <unordered_set>

/**
 * Funkcja pomocnicza do weryfikacji, czy dany rodzic jest poprawną
 * permutacją (albo w przybliżeniu poprawną) w zakresie 0..(size-1).
 * Jeżeli napotka out-of-range, to może:
 *  - wyzerować (uszkodzić) osobnika,
 *  - albo skorygować (np. clamp).
 * Tutaj dla przykładu – clamp.
 */
// Pomocnicza funkcja do sprawdzenia, czy chromosom jest prawidłowy (np. permutacja).
static bool validateAndClampChromosome(std::shared_ptr<std::vector<int>> indiv, int size) {
    if (!indiv) return false;
    if ((int)indiv->size() != size) {
        return false;
    }
    // W razie potrzeby można sprawdzić, czy wartości są unikatowe.
    for (int &gene : *indiv) {
        if (gene < 0 || gene >= size) {
            // "clamp" do zakresu
            if (gene < 0) gene = 0;
            if (gene >= size) gene = size - 1;
        }
    }
    return true;
}

// ========================================
// =            ONE POINT                =
// ========================================
std::vector<std::shared_ptr<std::vector<int>>> 
OnePointCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());

    // Zakładamy, że rozmiar = instance.getSpectrum().size()
    int requiredSize = (int)instance.getSpectrum().size();

    // Krzyżujemy parami: (0,1), (2,3), ...
    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i + 1];
        if (!p1 || !p2) {
            continue;
        }
        // weryfikacja
        bool ok1 = validateAndClampChromosome(p1, requiredSize);
        bool ok2 = validateAndClampChromosome(p2, requiredSize);
        if (!ok1 || !ok2) {
            // Skopiuj rodziców bez krzyżowania
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = (int)p1->size();

        // Losujemy punkt podziału
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, size - 1);
        int crossPoint = dist(gen);

        auto child1 = std::make_shared<std::vector<int>>(*p1);
        auto child2 = std::make_shared<std::vector<int>>(*p2);

        for (int idx = crossPoint; idx < size; idx++) {
            (*child1)[idx] = (*p2)[idx];
            (*child2)[idx] = (*p1)[idx];
        }

        offspring.push_back(child1);
        offspring.push_back(child2);
    }

    return offspring;
}

// ========================================
// =           ORDER CROSSOVER (OX)      =
// ========================================
std::vector<std::shared_ptr<std::vector<int>>> 
OrderCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());

    int requiredSize = (int)instance.getSpectrum().size();

    // co 2 rodziców
    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i+1];
        if (!p1 || !p2) {
            continue;
        }
        bool ok1 = validateAndClampChromosome(p1, requiredSize);
        bool ok2 = validateAndClampChromosome(p2, requiredSize);
        if (!ok1 || !ok2) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = (int)p1->size();

        // Losowanie segmentu
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<int> distIndex(0, size - 1);
        int start = distIndex(gen);
        int end   = distIndex(gen);
        if (start > end) std::swap(start, end);

        auto c1 = std::make_shared<std::vector<int>>(size, -1);
        auto c2 = std::make_shared<std::vector<int>>(size, -1);

        std::vector<bool> used1(size, false);
        std::vector<bool> used2(size, false);

        // Kopiujemy segment start..end
        for (int j = start; j <= end; j++) {
            (*c1)[j] = (*p1)[j];
            used1[(*p1)[j]] = true;

            (*c2)[j] = (*p2)[j];
            used2[(*p2)[j]] = true;
        }

        // wypełniamy resztę c1
        int curr1 = (end + 1) % size;
        int pos2  = (end + 1) % size;
        while (curr1 != start) {
            int val = (*p2)[pos2];
            if (!used1[val]) {
                (*c1)[curr1] = val;
                used1[val] = true;
                curr1 = (curr1 + 1) % size;
            }
            pos2 = (pos2 + 1) % size;
        }

        // wypełniamy resztę c2
        int curr2 = (end + 1) % size;
        int pos1  = (end + 1) % size;
        while (curr2 != start) {
            int val = (*p1)[pos1];
            if (!used2[val]) {
                (*c2)[curr2] = val;
                used2[val] = true;
                curr2 = (curr2 + 1) % size;
            }
            pos1 = (pos1 + 1) % size;
        }

        offspring.push_back(c1);
        offspring.push_back(c2);
    }
    return offspring;
}

// ========================================
// =       EDGE RECOMBINATION (ERX)      =
// ========================================
std::vector<std::shared_ptr<std::vector<int>>> 
EdgeRecombination::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());

    int requiredSize = (int)instance.getSpectrum().size();

    // parami
    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i+1];
        if (!p1 || !p2) {
            continue;
        }

        bool ok1 = validateAndClampChromosome(p1, requiredSize);
        bool ok2 = validateAndClampChromosome(p2, requiredSize);
        if (!ok1 || !ok2) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = (int)p1->size();
        std::vector<std::vector<int>> edgeMap(size);

        // Funkcja do budowy edge map
        auto addEdges = [&](std::shared_ptr<std::vector<int>> parent) {
            for (int idx = 0; idx < size; idx++) {
                int curr = (*parent)[idx];
                int next = (*parent)[(idx + 1) % size];
                int prev = (*parent)[(idx - 1 + size) % size];

                // proste sprawdzenie zakresu
                if (curr < 0 || curr >= size) continue;
                if (next < 0 || next >= size) continue;
                if (prev < 0 || prev >= size) continue;

                edgeMap[curr].push_back(prev);
                edgeMap[curr].push_back(next);
            }
        };

        addEdges(p1);
        addEdges(p2);

        // usuwamy duplikaty
        for (auto &edges : edgeMap) {
            std::sort(edges.begin(), edges.end());
            edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
        }

        // budowa childa
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<int> distStart(0, size - 1);

        // np. start od p1[0] lub losowo:
        int startGene = (*p1)[0]; 
        // int startGene = distStart(gen); // alternatywnie

        std::vector<int> child;
        child.reserve(size);
        child.push_back(startGene);

        std::vector<bool> used(size, false);
        used[startGene] = true;
        int current = startGene;

        while ((int)child.size() < size) {
            // Weźmy sąsiadów current
            const auto &neighbors = edgeMap[current];
            // odfiltrujmy nieużytych
            std::vector<int> validNbs;
            for (int nb : neighbors) {
                if (!used[nb]) {
                    validNbs.push_back(nb);
                }
            }

            if (validNbs.empty()) {
                // jeśli brak sąsiadów, wybieramy jakiś nieużyty gen z populacji
                std::vector<int> unusedGenes;
                for (int g = 0; g < size; g++) {
                    if (!used[g]) {
                        unusedGenes.push_back(g);
                    }
                }
                if (unusedGenes.empty()) {
                    // gotowe
                    break;
                }
                std::uniform_int_distribution<int> distU(0, (int)unusedGenes.size() - 1);
                int pick = unusedGenes[distU(gen)];
                child.push_back(pick);
                used[pick] = true;
                current = pick;
            } else {
                // losowo z validNbs
                std::uniform_int_distribution<int> distN(0, (int)validNbs.size() - 1);
                int pick = validNbs[distN(gen)];
                child.push_back(pick);
                used[pick] = true;
                current = pick;
            }

            // (opcjonalnie) dajmy break awaryjny jeśli cokolwiek idzie źle
            // np. if (child.size() > 2*size) break; 
            // – choć normalnie nie powinno się zdarzyć
        }

        // Mamy jednego potomka. W ERX często robi się 1 lub 2 potomków.
        // Stwórzmy jeszcze drugiego – np. odwrotnego:
        auto childPtr1 = std::make_shared<std::vector<int>>(child);
        auto revChild = child;
        std::reverse(revChild.begin(), revChild.end());
        auto childPtr2 = std::make_shared<std::vector<int>>(revChild);

        offspring.push_back(childPtr1);
        offspring.push_back(childPtr2);
    }

    return offspring;
}

// ========================================
// =        PARTIALLY MAPPED (PMX)        =
// ========================================
std::vector<std::shared_ptr<std::vector<int>>> 
PMXCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());

    int requiredSize = (int)instance.getSpectrum().size();

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i+1];
        if (!p1 || !p2) {
            continue;
        }
        bool ok1 = validateAndClampChromosome(p1, requiredSize);
        bool ok2 = validateAndClampChromosome(p2, requiredSize);
        if (!ok1 || !ok2) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = (int)p1->size();

        auto c1 = std::make_shared<std::vector<int>>(*p1);
        auto c2 = std::make_shared<std::vector<int>>(*p2);

        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, size - 1);
        int start = dist(gen);
        int end   = dist(gen);
        if (start > end) std::swap(start, end);

        // PMX mapowanie
        for (int pos = start; pos <= end; pos++) {
            int val1 = (*c1)[pos];
            int val2 = (*c2)[pos];
            // Szukamy pozycji val1 w c2
            // i pozycji val2 w c1
            // i zamieniamy
            for (int j = 0; j < size; j++) {
                if ((*c2)[j] == val1) {
                    (*c2)[j] = val2;
                    break;
                }
            }
            for (int j = 0; j < size; j++) {
                if ((*c1)[j] == val2) {
                    (*c1)[j] = val1;
                    break;
                }
            }
        }

        offspring.push_back(c1);
        offspring.push_back(c2);
    }

    return offspring;
}