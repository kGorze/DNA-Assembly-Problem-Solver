//
// Created by konrad_guest on 09/01/2025.
// SMART
#include "metaheuristics/population_cache.h"
#include "metaheuristics/fitness.h"
#include <algorithm>
#include <chrono>
#include <mutex>

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
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation
) {
    std::unique_lock<std::mutex> lock(m_mutex);
    
    // Try to find in cache
    auto it = m_cache.find(solution);
    if (it != m_cache.end()) {
        return it->second;
    }
    
    // Calculate and cache
    double value = fitness->calculateFitness(solution, instance, representation);
    m_cache[solution] = value;
    return value;
}

void CachedPopulation::updatePopulation(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation
) {
    std::unique_lock<std::mutex> lock(m_mutex);
    
    // Clear old entries
    m_cache.clear();
    
    // Cache new population
    for (const auto& solution : population) {
        m_cache[solution] = fitness->calculateFitness(solution, instance, representation);
    }
}

void CachedPopulation::clear() {
    std::unique_lock<std::mutex> lock(m_mutex);
    m_cache.clear();
}

size_t CachedPopulation::size() const {
    return cache.size();
}
