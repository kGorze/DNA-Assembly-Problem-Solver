//
// Created by konrad_guest on 09/01/2025.
// SMART
#include "../../include/metaheuristics/population_cache_impl.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <chrono>
#include <mutex>
#include <sstream>
#include <limits>

double SimplePopulationCache::getOrCalculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance)
{
    if (!individual) {
        LOG_ERROR("Cannot calculate fitness for null individual");
        return -std::numeric_limits<double>::infinity();
    }
    
    std::lock_guard<std::mutex> lock(cacheMutex);
    
    if (!m_fitness) {
        LOG_ERROR("Fitness calculator not initialized");
        return -std::numeric_limits<double>::infinity();
    }
    
    try {
        std::string key = individual->toString();
        
        // Check cache first
        auto it = cache.find(key);
        if (it != cache.end()) {
            return it->second;
        }
        
        // Calculate fitness
        double fitness = m_fitness->evaluate(individual, instance);
        
        // Cache result
        if (cache.size() >= maxCacheSize) {
            cleanupCache();
        }
        cache[key] = fitness;
        
        return fitness;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to calculate fitness: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

void SimplePopulationCache::updatePopulation(
    const std::vector<std::shared_ptr<Individual>>& population)
{
    if (population.empty()) {
        return;
    }
    
    std::lock_guard<std::mutex> lock(cacheMutex);
    
    // Clear old entries if cache is too large
    if (cache.size() + population.size() > maxCacheSize) {
        cleanupCache();
    }
    
    // Update cache with new population
    for (const auto& individual : population) {
        if (individual) {
            try {
                std::string key = individual->toString();
                if (auto it = cache.find(key); it == cache.end()) {
                    cache[key] = individual->getFitness();
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Failed to update cache for individual: " + std::string(e.what()));
                continue;
            }
        }
    }
}

void SimplePopulationCache::clear() {
    std::lock_guard<std::mutex> lock(cacheMutex);
    cache.clear();
}

void SimplePopulationCache::cleanupCache() {
    if (cache.size() <= maxCacheSize / 2) return;
    
    try {
        std::vector<std::pair<std::string, double>> entries(
            cache.begin(), cache.end());
        
        // Keep only the most recent half
        size_t keepCount = maxCacheSize / 2;
        cache.clear();
        for (size_t i = entries.size() - keepCount; i < entries.size(); ++i) {
            cache.insert(entries[i]);
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to cleanup cache: " + std::string(e.what()));
        // If cleanup fails, just clear the cache entirely
        cache.clear();
    }
}
