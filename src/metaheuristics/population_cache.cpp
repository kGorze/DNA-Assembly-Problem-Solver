//
// Created by konrad_guest on 09/01/2025.
// SMART
#include "../include/metaheuristics/population_cache_impl.h"
#include "../include/utils/logging.h"
#include <algorithm>
#include <chrono>
#include <mutex>
#include <sstream>

double SimplePopulationCache::getOrCalculateFitness(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<const IFitness> fitness,
    std::shared_ptr<IRepresentation> representation
) {
    if (!solution || solution->empty()) {
        LOG_WARNING("Attempted to calculate fitness for null or empty solution");
        return 0.0;
    }

    std::stringstream ss;
    for (int gene : *solution) {
        ss << gene << ",";
    }
    std::string key = ss.str();

    {
        std::shared_lock<std::shared_mutex> readLock(m_mutex);
        auto it = m_cache.find(key);
        if (it != m_cache.end()) {
            return it->second;
        }
    }

    double calculatedFitness = fitness->calculateFitness(solution, instance, representation);

    {
        std::unique_lock<std::shared_mutex> writeLock(m_mutex);
        if (m_cache.size() >= MAX_CACHE_SIZE) {
            m_cache.clear();
            LOG_INFO("Cache size limit reached, clearing cache");
        }
        m_cache[key] = calculatedFitness;
    }

    return calculatedFitness;
}

void SimplePopulationCache::updatePopulation(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const DNAInstance& instance,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IRepresentation> representation
) {
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_cache.clear();
    
    for (const auto& individual : population) {
        if (individual) {
            std::stringstream ss;
            for (int gene : *individual) {
                ss << gene << ",";
            }
            m_cache[ss.str()] = fitness->calculateFitness(individual, instance, representation);
        }
    }
}

void SimplePopulationCache::clear() {
    std::unique_lock<std::shared_mutex> lock(m_mutex);
    m_cache.clear();
    LOG_INFO("Population cache cleared");
}
