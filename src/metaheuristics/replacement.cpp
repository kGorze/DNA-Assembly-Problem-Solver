//
// Created by konrad_guest on 07/01/2025.
//

#include "metaheuristics/replacement.h"
#include <algorithm>

// =========== FullReplacement ===========
std::vector<void*> 
FullReplacement::replace(const std::vector<void*>        &oldPop,
                         const std::vector<void*>        &offspring,
                         const DNAInstance               &instance,
                         std::shared_ptr<IFitness>        fitness,
                         std::shared_ptr<IRepresentation> representation)
{
    // trivial: discard oldPop, return offspring
    return offspring;
}

// =========== PartialReplacement ===========
class PartialReplacement::PopulationCache {
public:
    struct Individual {
        void* ptr;
        double fitness;
        
        Individual(void* p = nullptr, double f = -std::numeric_limits<double>::infinity()) 
            : ptr(p), fitness(f) {}
        
        bool operator<(const Individual& other) const {
            return fitness > other.fitness; // Sort in descending order
        }
    };
    
    std::vector<Individual> sorted;
    bool needsSort = true;
    
    void clear() {
        sorted.clear();
        needsSort = true;
    }
};

PartialReplacement::PartialReplacement(double replacementRatio) 
    : m_replacementRatio(replacementRatio)
    , m_cache(std::make_unique<PopulationCache>()) {}

PartialReplacement::~PartialReplacement() = default;

std::vector<void*> 
PartialReplacement::replace(const std::vector<void*>        &oldPop,
                            const std::vector<void*>        &offspring,
                            const DNAInstance               &instance,
                            std::shared_ptr<IFitness>        fitness,
                            std::shared_ptr<IRepresentation> representation)
{
    // Ensure we have valid input populations
    if (oldPop.empty() && offspring.empty()) {
        return std::vector<void*>();
    }

    const size_t oldSize = oldPop.size();
    const size_t offSize = offspring.size();
    const size_t targetSize = std::max(oldSize, offSize); // Maintain population size
    
    // Calculate how many parents to keep
    const size_t desiredParents = static_cast<size_t>(targetSize * (1.0 - m_replacementRatio));
    
    std::vector<void*> newPop;
    newPop.reserve(targetSize);
    
    // First, evaluate and sort old population
    std::vector<PopulationCache::Individual> allIndividuals;
    allIndividuals.reserve(oldSize + offSize);
    
    // Add old population individuals
    for (auto* ind : oldPop) {
        if (ind != nullptr) {
            double fit = fitness->evaluate(ind, instance, representation);
            allIndividuals.emplace_back(ind, fit);
        }
    }
    
    // Add offspring individuals
    for (auto* ind : offspring) {
        if (ind != nullptr) {
            double fit = fitness->evaluate(ind, instance, representation);
            allIndividuals.emplace_back(ind, fit);
        }
    }
    
    // Sort all individuals by fitness
    std::sort(allIndividuals.begin(), allIndividuals.end());
    
    // Take the best individuals up to target size
    size_t individualsToTake = std::min(targetSize, allIndividuals.size());
    
    for (size_t i = 0; i < individualsToTake; ++i) {
        auto* orig = static_cast<std::vector<int>*>(allIndividuals[i].ptr);
        if (orig != nullptr) {
            auto* copy = new std::vector<int>(*orig);
            newPop.push_back(copy);
        }
    }
    
    // Clean up old individuals that weren't selected
    for (size_t i = individualsToTake; i < allIndividuals.size(); ++i) {
        auto* ptr = static_cast<std::vector<int>*>(allIndividuals[i].ptr);
        delete ptr;
    }
    
    return newPop;
}
