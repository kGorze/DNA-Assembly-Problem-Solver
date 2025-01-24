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
// =            ONE POINT                =
// ========================================
std::vector<std::shared_ptr<std::vector<int>>> 
OnePointCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    LOG_DEBUG("Starting crossover operation with " + std::to_string(parents.size()) + " parents");
    
    if (parents.size() < 2) {
        LOG_ERROR("Insufficient parents for crossover. Need at least 2, got " + std::to_string(parents.size()));
        return parents;
    }
    
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());
    
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    for (size_t i = 0; i < parents.size(); i += 2) {
        if (i + 1 >= parents.size()) {
            LOG_WARNING("Odd number of parents, copying last parent directly");
            offspring.push_back(std::make_shared<std::vector<int>>(*parents[i]));
            continue;
        }
        
        if (dist(rng) < m_crossoverRate) {
            DEBUG_LOG("Performing crossover between parents " + std::to_string(i) + 
                     " and " + std::to_string(i + 1));
            
            auto children = performCrossover(parents[i], parents[i + 1], instance);
            
            if (children.first && children.second) {
                offspring.push_back(children.first);
                offspring.push_back(children.second);
            } else {
                LOG_ERROR("Crossover failed, using parents instead");
                offspring.push_back(std::make_shared<std::vector<int>>(*parents[i]));
                offspring.push_back(std::make_shared<std::vector<int>>(*parents[i + 1]));
            }
        } else {
            DEBUG_LOG("Skipping crossover for parents " + std::to_string(i) + 
                     " and " + std::to_string(i + 1));
            offspring.push_back(std::make_shared<std::vector<int>>(*parents[i]));
            offspring.push_back(std::make_shared<std::vector<int>>(*parents[i + 1]));
        }
    }
    
    LOG_INFO("Crossover complete. Generated " + std::to_string(offspring.size()) + " offspring");
    return offspring;
}

std::pair<std::shared_ptr<std::vector<int>>, std::shared_ptr<std::vector<int>>> 
OnePointCrossover::performCrossover(const std::shared_ptr<std::vector<int>>& parent1,
                                  const std::shared_ptr<std::vector<int>>& parent2,
                                  const DNAInstance& instance) {
    auto child1 = std::make_shared<std::vector<int>>(*parent1);
    auto child2 = std::make_shared<std::vector<int>>(*parent2);
    
    std::random_device rd;
    std::mt19937 gen(rd());
    int crossPoint = std::uniform_int_distribution<int>(0, parent1->size()-1)(gen);
    
    for (size_t i = crossPoint; i < parent1->size(); ++i) {
        (*child1)[i] = (*parent2)[i];
        (*child2)[i] = (*parent1)[i];
    }
    
    return {child1, child2};
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
    const size_t spectrumSize = instance.getSpectrum().size();
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    
    if (parents.size() < 2) {
        std::cerr << "[ERROR OrderCrossover] Not enough parents for crossover\n";
        return offspring;
    }

    for (size_t i = 0; i < parents.size() - 1; i += 2) {
        if (!parents[i] || !parents[i + 1] || 
            parents[i]->size() != spectrumSize || 
            parents[i + 1]->size() != spectrumSize) {
            std::cerr << "[ERROR OrderCrossover] Invalid parent sizes or null parents\n";
            continue;
        }

        auto child = performOrderCrossover(*parents[i], *parents[i + 1], spectrumSize);
        if (isValidPermutation(child, spectrumSize)) {
            offspring.push_back(std::make_shared<std::vector<int>>(child));
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

    int requiredSize = parents.empty() ? 0 : (int)parents[0]->size();

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
            for (auto& child : offspring) {
                if (child) {
                    for (int gene : *child) {
                        if (gene < 0 || gene >= requiredSize) {
                            std::cerr << "[ERROR] child has out-of-range gene = " << gene 
                                      << " (valid range: 0-" << (requiredSize-1) << ")" << std::endl;
                        }
                    }
                }
            }
        } else {
            // Fallback do kopii rodziców
            std::cerr << "[ERROR] Invalid offspring from EdgeRecombination!" << std::endl;
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

    int requiredSize = parents.empty() ? 0 : (int)parents[0]->size();

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i+1];
        // **Debug Log**: Sprawdzenie rozmiarów rodziców
        std::cout << "[DEBUG PMX] Parent1 size: " << p1->size() 
                  << ", Parent2 size: " << p2->size() << "\n";

        if (!p1 || !p2) continue;

        if (!isValidPermutation(*p1, requiredSize) || !isValidPermutation(*p2, requiredSize)) {
            std::cerr << "[ERROR PMX] Invalid permutation size for parents\n";

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
            // **Debug Log**: Sprawdzenie potomków
            for(auto &child : {c1, c2}) {
                std::cout << "[DEBUG PMX] Child size: " << child->size() << ", genes: ";
                for(auto gene : *child){
                    std::cout << gene << " ";
                    if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                        std::cerr << "\n[ERROR PMX] Child has out-of-range gene: " << gene << "\n";
                    }
                }
                std::cout << "\n";
            }
        } else {
            std::cerr << "[ERROR] Invalid offspring from PMXCrossover!" << std::endl;
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
        }
    }

    return offspring;
}

DistancePreservingCrossover::DistanceMatrix::DistanceMatrix(const std::vector<int>& perm) {
    int size = perm.size();
    // Zamiast przechowywać pełną macierz, przechowujemy tylko pozycje elementów
    distances.resize(size);
    for (int i = 0; i < size; i++) {
        distances[perm[i]] = i;
    }
}

int DistancePreservingCrossover::DistanceMatrix::getDistance(int from, int to) const {
    int fromPos = distances[from];
    int toPos = distances[to];
    int size = distances.size();
    
    // Oblicz najkrótszą odległość w cyklu
    int dist = (toPos - fromPos + size) % size;
    int revDist = (fromPos - toPos + size) % size;
    return std::min(dist, revDist);
}

std::vector<std::shared_ptr<std::vector<int>>> 
DistancePreservingCrossover::crossover(
    const std::vector<std::shared_ptr<std::vector<int>>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    offspring.reserve(parents.size());

    int requiredSize = parents.empty() ? 0 : (int)parents[0]->size();

    std::mt19937 gen(std::random_device{}());
    std::vector<bool> used(requiredSize);
    std::vector<std::pair<int, double>> candidates;
    candidates.reserve(requiredSize);

    for (size_t i = 0; i + 1 < parents.size(); i += 2) {
        auto p1 = parents[i];
        auto p2 = parents[i + 1];
        if (!p1 || !p2) continue;

        if (!isValidPermutation(*p1, requiredSize) || !isValidPermutation(*p2, requiredSize)) {
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
            continue;
        }

        // Tworzymy macierze odległości dla obu rodziców - teraz O(n) zamiast O(n^2)
        DistanceMatrix dm1(*p1);
        DistanceMatrix dm2(*p2);

        // Tworzymy potomka
        auto child = std::make_shared<std::vector<int>>(requiredSize);
        std::fill(used.begin(), used.end(), false);

        // Wybierz pierwszy element
        (*child)[0] = (*p1)[0];
        used[(*child)[0]] = true;

        // Zoptymalizowana pętla główna
        for (int pos = 1; pos < requiredSize; pos++) {
            double bestScore = std::numeric_limits<double>::infinity();
            int bestCandidate = -1;
            
            // Zamiast sprawdzać wszystkie możliwe kombinacje, weźmiemy próbkę
            int sampleSize = std::min(10, requiredSize - pos); // Sprawdzamy max 10 kandydatów
            std::vector<int> unusedIndices;
            for (int j = 0; j < requiredSize; j++) {
                if (!used[j]) {
                    unusedIndices.push_back(j);
                }
            }
            std::shuffle(unusedIndices.begin(), unusedIndices.end(), gen);
            
            // Sprawdź tylko wybraną próbkę kandydatów
            for (int idx = 0; idx < sampleSize && idx < (int)unusedIndices.size(); idx++) {
                int j = unusedIndices[idx];
                double score = 0.0;
                
                // Sprawdź tylko względem kilku ostatnich umieszczonych elementów
                int checkBack = std::min(5, pos); // Sprawdzamy max 5 poprzednich elementów
                for (int k = pos - checkBack; k < pos; k++) {
                    int placedElement = (*child)[k];
                    score += std::abs(dm1.getDistance(placedElement, j) - (pos - k)) +
                            std::abs(dm2.getDistance(placedElement, j) - (pos - k));
                }
                
                if (score < bestScore) {
                    bestScore = score;
                    bestCandidate = j;
                }
            }
            
            (*child)[pos] = bestCandidate;
            used[bestCandidate] = true;
        }

        if (isValidPermutation(*child, requiredSize)) {
            offspring.push_back(child);
            
            auto child2 = std::make_shared<std::vector<int>>(*child);
            std::reverse(child2->begin(), child2->end());
            offspring.push_back(child2);
            for (auto& child : offspring) {
                if (child) {
                    for (int gene : *child) {
                        if (gene < 0 || gene >= requiredSize) {
                            std::cerr << "[ERROR] child has out-of-range gene = " << gene 
                                      << " (valid range: 0-" << (requiredSize-1) << ")" << std::endl;
                        }
                    }
                }
            }
        } else {
            std::cerr << "[ERROR] Invalid offspring from DistancePreservingCrossover!" << std::endl;
            offspring.push_back(std::make_shared<std::vector<int>>(*p1));
            offspring.push_back(std::make_shared<std::vector<int>>(*p2));
        }
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
    
    // Wybierz losowy segment
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, size - 1);
    int start = dist(gen);
    int end = dist(gen);
    if (start > end) std::swap(start, end);
    
    // Skopiuj segment z pierwszego rodzica
    for (int i = start; i <= end; i++) {
        child[i] = parent1[i];
        used[parent1[i]] = true;
    }
    
    // Wypełnij pozostałe pozycje genami z drugiego rodzica
    int j = (end + 1) % size;
    for (size_t i = 0; i < size; i++) {
        int pos = (end + 1 + i) % size;
        if (child[pos] == -1) {
            // Znajdź następny nieużyty gen z parent2
            while (used[parent2[j]]) {
                j = (j + 1) % size;
            }
            child[pos] = parent2[j];
            used[parent2[j]] = true;
            j = (j + 1) % size;
        }
    }
    
    return child;
}