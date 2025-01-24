//
// Created by konrad_guest on 09/01/2025.
// SMART
#include "metaheuristics/population_cache.h"
#include "metaheuristics/fitness.h"
#include <algorithm>
#include <chrono>

uint64_t CachedPopulation::computeHash(std::shared_ptr<std::vector<int>> individual) {
    return ZobristHasher::getInstance().hashPermutation(*individual);
}

void CachedPopulation::evictOldEntries() {
    if (cache.size() <= maxCacheSize) return;
    
    std::vector<std::pair<uint64_t, std::chrono::steady_clock::time_point>> entries;
    entries.reserve(cache.size());
    
    for (const auto& entry : cache) {
        entries.emplace_back(entry.first, entry.second.lastAccess);
    }
    
    std::sort(entries.begin(), entries.end(),
             [](const auto& a, const auto& b) {
                 return a.second < b.second;
             });
    
    size_t toRemove = cache.size() / 5; // remove 20% oldest
    for (size_t i = 0; i < toRemove; ++i) {
        cache.erase(entries[i].first);
    }
}

double CachedPopulation::getOrCalculateFitness(
    std::shared_ptr<std::vector<int>> individual,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation) 
{
    uint64_t hash = computeHash(individual);
    
    auto it = cache.find(hash);
    if (it != cache.end()) {
        it->second.lastAccess = std::chrono::steady_clock::now();
        return it->second.fitness;
    }
    
    double fitnessValue = fitness->evaluate(individual, instance, representation);
    
    cache[hash] = {
        fitnessValue,
        hash,
        std::chrono::steady_clock::now()
    };
    
    evictOldEntries();
    
    return fitnessValue;
}

void CachedPopulation::updatePopulation(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation) 
{
    for (auto& individual : population) {
        getOrCalculateFitness(individual, instance, fitness, representation);
    }
}

void CachedPopulation::clear() {
    cache.clear();
}

size_t CachedPopulation::size() const {
    return cache.size();
}
