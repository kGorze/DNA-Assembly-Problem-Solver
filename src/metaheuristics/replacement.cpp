//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/replacement.h"
#include <algorithm>

// =========== FullReplacement ===========

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

    struct Individual {
        std::shared_ptr<std::vector<int>> ptr;
        double fitnessVal;
        
        Individual(std::shared_ptr<std::vector<int>> p = nullptr, double f = -std::numeric_limits<double>::infinity()) 
            : ptr(p), fitnessVal(f) {}
        
        bool operator<(const Individual& other) const {
            return fitnessVal > other.fitnessVal; // descending order
        }
    };
    
    std::vector<Individual> allIndividuals;
    allIndividuals.reserve(oldSize + offSize);
    
    for (auto &ind : oldPop) {
        if (ind) {
            double fit = m_fitnessCache->getOrCalculateFitness(ind, instance, fitness, representation);
            allIndividuals.emplace_back(ind, fit);
        }
    }
    
    for (auto &ind : offspring) {
        if (ind) {
            double fit = m_fitnessCache->getOrCalculateFitness(ind, instance, fitness, representation);
            allIndividuals.emplace_back(ind, fit);
        }
    }
    
    std::sort(allIndividuals.begin(), allIndividuals.end());
    
    std::vector<std::shared_ptr<std::vector<int>>> newPop;
    newPop.reserve(targetSize);
    size_t individualsToTake = std::min(targetSize, allIndividuals.size());
    
    for (size_t i = 0; i < individualsToTake; ++i) {
        newPop.push_back(allIndividuals[i].ptr);
    }

    return newPop;
}
