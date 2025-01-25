#pragma once

#include "metaheuristics/individual.h"
#include "dna/dna_instance.h"
#include <memory>
#include <unordered_map>
#include <mutex>

class IPopulationCache {
public:
    virtual ~IPopulationCache() = default;
    virtual double getOrCalculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) = 0;
    virtual void clear() = 0;
};

class SimplePopulationCache : public IPopulationCache {
public:
    static constexpr size_t MAX_CACHE_SIZE = 10000;

    SimplePopulationCache() : m_maxCacheSize(MAX_CACHE_SIZE) {}
    explicit SimplePopulationCache(size_t cacheSize) 
        : m_maxCacheSize(std::max(size_t(1), cacheSize)) {}

    double getOrCalculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) override {
        if (!individual) return 0.0;

        std::lock_guard<std::mutex> lock(m_mutex);
        
        // Try to find in cache
        auto it = m_cache.find(individual->toString());
        if (it != m_cache.end()) {
            return it->second;
        }

        // Calculate fitness
        double fitness = individual->getFitness();

        // Add to cache if not full
        if (m_cache.size() < m_maxCacheSize) {
            m_cache[individual->toString()] = fitness;
        }

        return fitness;
    }

    void clear() override {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_cache.clear();
    }

private:
    const size_t m_maxCacheSize;
    std::unordered_map<std::string, double> m_cache;
    mutable std::mutex m_mutex;
}; 