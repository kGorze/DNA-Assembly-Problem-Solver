#pragma once

#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include <unordered_map>
#include <mutex>
#include <memory>
#include <sstream>
#include <string>
#include <limits>

class SimplePopulationCache : public IPopulationCache {
public:
    static constexpr size_t MAX_CACHE_SIZE = 10000;

    SimplePopulationCache() : maxCacheSize(MAX_CACHE_SIZE), m_fitness(nullptr) {}
    
    explicit SimplePopulationCache(size_t cacheSize) 
        : maxCacheSize(cacheSize), m_fitness(nullptr) {}

    double getOrCalculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance) override;

    void updatePopulation(
        const std::vector<std::shared_ptr<Individual>>& population) override;

    void clear() override;

    size_t size() const {
        std::lock_guard<std::mutex> lock(cacheMutex);
        return cache.size();
    }

    void setFitnessCalculator(std::shared_ptr<IFitness> fitness) {
        std::lock_guard<std::mutex> lock(cacheMutex);
        m_fitness = std::move(fitness);
    }

private:
    std::unordered_map<std::string, double> cache;
    mutable std::mutex cacheMutex;
    const size_t maxCacheSize;
    std::shared_ptr<IFitness> m_fitness;
    
    void cleanupCache();
}; 