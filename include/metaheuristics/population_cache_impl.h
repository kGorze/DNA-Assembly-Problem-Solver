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

    SimplePopulationCache() : m_maxCacheSize(MAX_CACHE_SIZE), m_fitness(nullptr) {}
    
    explicit SimplePopulationCache(size_t cacheSize) 
        : m_maxCacheSize(std::max(size_t(1), cacheSize)), m_fitness(nullptr) {}

    double getOrCalculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance) override;

    void updatePopulation(
        const std::vector<std::shared_ptr<Individual>>& population) override;

    void clear() override {
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        m_cache.clear();
    }

    void reserve(size_t size) override {
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        m_cache.reserve(size);
    }

    void add(const std::shared_ptr<Individual>& individual) override {
        if (!individual) return;
        
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        try {
            std::string key = individual->toString();
            if (m_cache.size() >= m_maxCacheSize) {
                cleanupCache();
            }
            m_cache[key] = individual->getFitness();
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to add individual to cache: " + std::string(e.what()));
        }
    }

    bool contains(const std::shared_ptr<Individual>& individual) const override {
        if (!individual) return false;
        
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        try {
            return m_cache.find(individual->toString()) != m_cache.end();
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to check cache: " + std::string(e.what()));
            return false;
        }
    }

    size_t size() const {
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        return m_cache.size();
    }

    void setFitnessCalculator(std::shared_ptr<IFitness> fitness) {
        std::lock_guard<std::mutex> lock(m_cacheMutex);
        m_fitness = std::move(fitness);
    }

private:
    std::unordered_map<std::string, double> m_cache;
    mutable std::mutex m_cacheMutex;
    const size_t m_maxCacheSize;
    std::shared_ptr<IFitness> m_fitness;
    
    void cleanupCache();
}; 