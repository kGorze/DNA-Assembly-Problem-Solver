//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/replacement.h"
#include <algorithm>
#include <iostream>

std::vector<std::shared_ptr<std::vector<int>>>
FullReplacement::replace(const std::vector<std::shared_ptr<std::vector<int>>> &oldPop,
                         const std::vector<std::shared_ptr<std::vector<int>>> &offspring,
                         const DNAInstance &instance,
                         std::shared_ptr<IFitness> fitness,
                         std::shared_ptr<IRepresentation> representation)
{
    // trivial: discard oldPop, return offspring
    return offspring;
}

// =========== PartialReplacement ===========

PartialReplacement::PartialReplacement(double replacementRatio, std::shared_ptr<IPopulationCache> cache) 
    : m_replacementRatio(replacementRatio)
    , m_fitnessCache(cache) {}

PartialReplacement::~PartialReplacement() = default;

std::vector<std::shared_ptr<std::vector<int>>> 
PartialReplacement::replace(const std::vector<std::shared_ptr<std::vector<int>>> &oldPop,
                            const std::vector<std::shared_ptr<std::vector<int>>> &offspring,
                            const DNAInstance &instance,
                            std::shared_ptr<IFitness> fitness,
                            std::shared_ptr<IRepresentation> representation)
{
    if (oldPop.empty() && offspring.empty()) {
        return {};
    }

    const size_t oldSize = oldPop.size();
    const size_t offSize = offspring.size();
    const size_t targetSize = std::max(oldSize, offSize);
    
    // Calculate required length based on n and k (same as in PermutationRepresentation)
    int k = instance.getK();
    int n = instance.getN();
    int requiredLength = n - k + 1;  // This is the correct length for Permutation representation

    struct Individual {
        std::shared_ptr<std::vector<int>> ptr;
        double fitnessVal;
        
        Individual(std::shared_ptr<std::vector<int>> p = nullptr, double f = -1e9) 
            : ptr(p), fitnessVal(f) {}
    };
    
    std::vector<Individual> allIndividuals;
    allIndividuals.reserve(oldSize + offSize);
    
    // Lambda for size verification
    auto isValidSize = [&](std::shared_ptr<std::vector<int>> ind){
        return (ind && (int)ind->size() == requiredLength);
    };

    // 1. Zbieramy stare + nowe
    // 2. Odfiltrowujemy osobniki z niewłaściwą długością, bo i tak spowodują błąd
    // 3. Obliczamy fitness
    auto addWithFitness = [&](const std::vector<std::shared_ptr<std::vector<int>>>& vec) {
        for (auto &ind : vec) {
            if (!ind) continue;
            if (!isValidSize(ind)) {
                // Możesz logować ostrzeżenie
                continue;
            }
            double fit = m_fitnessCache->getOrCalculateFitness(ind, instance, fitness, representation);
            allIndividuals.push_back(Individual(ind, fit));
        }
    };

    addWithFitness(oldPop);
    addWithFitness(offspring);
    
    // Jeśli wszystko odrzucone, zwracamy cokolwiek
    if (allIndividuals.empty()) {
        // lub stwórzmy nową populację
        return representation->initializePopulation(targetSize, instance);
    }
    
    // Sortujemy malejąco wg fitnessVal
    std::sort(allIndividuals.begin(), allIndividuals.end(),
        [](const Individual& a, const Individual& b){
            return a.fitnessVal > b.fitnessVal;
        }
    );
    
    // Zostawiamy top N (N=targetSize)
    std::vector<std::shared_ptr<std::vector<int>>> newPop;
    newPop.reserve(targetSize);
    size_t take = std::min(targetSize, allIndividuals.size());
    for (size_t i = 0; i < take; i++) {
        newPop.push_back(allIndividuals[i].ptr);
    }

    return newPop;
}
