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

static bool isValidPermutation(const std::vector<int>& vec, int size) {
    if ((int)vec.size() != size) return false;
    std::vector<bool> used(size, false);
    for (int val : vec) {
        if (val < 0 || val >= size || used[val]) return false;
        used[val] = true;
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

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i+1];
        if (!p1 || !p2) continue;

        // Sprawdzamy czy rodzice są poprawnymi permutacjami
        if (!isValidPermutation(*p1, requiredSize) || !isValidPermutation(*p2, requiredSize)) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = requiredSize;
        
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<int> distIndex(0, size - 1);
        int start = distIndex(gen);
        int end = distIndex(gen);
        if (start > end) std::swap(start, end);

        // Inicjalizacja potomków
        auto c1 = std::make_shared<std::vector<int>>(size, -1);
        auto c2 = std::make_shared<std::vector<int>>(size, -1);

        std::vector<bool> used1(size, false);
        std::vector<bool> used2(size, false);

        // Kopiowanie segmentu
        for (int j = start; j <= end; j++) {
            (*c1)[j] = (*p1)[j];
            used1[(*p1)[j]] = true;

            (*c2)[j] = (*p2)[j];
            used2[(*p2)[j]] = true;
        }

        // Wypełnianie pozostałych pozycji
        int remainingCount1 = size - (end - start + 1);
        int remainingCount2 = remainingCount1;
        int curr1 = (end + 1) % size;
        int curr2 = (end + 1) % size;
        int pos1 = (end + 1) % size;
        int pos2 = (end + 1) % size;

        // Bezpieczne wypełnianie c1
        while (remainingCount1 > 0) {
            int val = (*p2)[pos2];
            if (!used1[val]) {
                (*c1)[curr1] = val;
                used1[val] = true;
                remainingCount1--;
                curr1 = (curr1 + 1) % size;
            }
            pos2 = (pos2 + 1) % size;
        }

        // Bezpieczne wypełnianie c2
        while (remainingCount2 > 0) {
            int val = (*p1)[pos1];
            if (!used2[val]) {
                (*c2)[curr2] = val;
                used2[val] = true;
                remainingCount2--;
                curr2 = (curr2 + 1) % size;
            }
            pos1 = (pos1 + 1) % size;
        }

        // Weryfikacja potomków
        if (isValidPermutation(*c1, size) && isValidPermutation(*c2, size)) {
            offspring.push_back(c1);
            offspring.push_back(c2);
        } else {
            // Jeśli coś poszło nie tak, zwracamy kopie rodziców
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
        }
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

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i+1];
        if (!p1 || !p2) continue;

        if (!isValidPermutation(*p1, requiredSize) || !isValidPermutation(*p2, requiredSize)) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = requiredSize;
        std::vector<std::unordered_set<int>> edgeMap(size);

        // Budowa edge map z użyciem set dla unikatowych wartości
        for (auto parent : {p1, p2}) {
            for (int idx = 0; idx < size; idx++) {
                int curr = (*parent)[idx];
                int next = (*parent)[(idx + 1) % size];
                int prev = (*parent)[(idx - 1 + size) % size];
                
                edgeMap[curr].insert(prev);
                edgeMap[curr].insert(next);
            }
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        
        // Tworzymy pierwszego potomka
        std::vector<bool> used(size, false);
        std::vector<int> child;
        child.reserve(size);
        
        // Zaczynamy od pierwszego genu pierwszego rodzica
        int current = (*p1)[0];
        child.push_back(current);
        used[current] = true;

        // Główna pętla budowania potomka
        while (child.size() < size) {
            std::vector<std::pair<int, int>> candidates; // para (gen, liczba nieużytych sąsiadów)
            
            // Sprawdzamy wszystkich nieużytych sąsiadów
            for (int neighbor : edgeMap[current]) {
                if (!used[neighbor]) {
                    int unusedNeighborCount = 0;
                    for (int nb : edgeMap[neighbor]) {
                        if (!used[nb]) unusedNeighborCount++;
                    }
                    candidates.push_back({neighbor, unusedNeighborCount});
                }
            }

            // Jeśli nie ma kandydatów, wybieramy spośród wszystkich nieużytych
            if (candidates.empty()) {
                for (int j = 0; j < size; j++) {
                    if (!used[j]) {
                        candidates.push_back({j, 0});
                    }
                }
            }

            // Wybieramy następny gen
            if (!candidates.empty()) {
                // Sortujemy po liczbie nieużytych sąsiadów (heurystyka)
                std::sort(candidates.begin(), candidates.end(),
                    [](const auto& a, const auto& b) { return a.second < b.second; });
                
                current = candidates[0].first;
                child.push_back(current);
                used[current] = true;
            } else {
                break; // Zabezpieczenie przed nieskończoną pętlą
            }
        }

        // Weryfikacja i dodanie potomków
        if (isValidPermutation(child, size)) {
            auto childPtr1 = std::make_shared<std::vector<int>>(child);
            // Drugi potomek jako odwrócona kopia pierwszego
            auto revChild = child;
            std::reverse(revChild.begin(), revChild.end());
            auto childPtr2 = std::make_shared<std::vector<int>>(revChild);
            
            offspring.push_back(childPtr1);
            offspring.push_back(childPtr2);
        } else {
            // Fallback do kopii rodziców
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
        }
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
        if (!p1 || !p2) continue;

        if (!isValidPermutation(*p1, requiredSize) || !isValidPermutation(*p2, requiredSize)) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        int size = requiredSize;

        // Tworzymy mapy pozycji dla obu rodziców
        std::vector<int> pos1(size), pos2(size);
        for (int j = 0; j < size; j++) {
            pos1[(*p1)[j]] = j;
            pos2[(*p2)[j]] = j;
        }

        auto c1 = std::make_shared<std::vector<int>>(*p1);
        auto c2 = std::make_shared<std::vector<int>>(*p2);

        // Losowanie punktów cięcia
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, size - 1);
        int start = dist(gen);
        int end = dist(gen);
        if (start > end) std::swap(start, end);

        // Wykonujemy mapowanie PMX
        for (int j = start; j <= end; j++) {
            int val1 = (*c1)[j];
            int val2 = (*c2)[j];
            
            std::swap((*c1)[j], (*c2)[j]);
            std::swap((*c1)[pos1[val2]], (*c2)[pos2[val1]]);
            
            // Aktualizujemy mapy pozycji
            int tempPos = pos1[val1];
            pos1[val1] = pos1[val2];
            pos1[val2] = tempPos;
            
            tempPos = pos2[val1];
            pos2[val1] = pos2[val2];
            pos2[val2] = tempPos;
        }

        // Weryfikacja potomków
        if (isValidPermutation(*c1, size) && isValidPermutation(*c2, size)) {
            offspring.push_back(c1);
            offspring.push_back(c2);
        } else {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
        }
    }

    return offspring;
}