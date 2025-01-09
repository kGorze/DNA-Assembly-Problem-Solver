//
// Created by konrad_guest on 09/01/2025.
//

#ifndef POPULATION_CACHE_H
#define POPULATION_CACHE_H

#include <vector>
#include <unordered_map>


#include "generator/dna_generator.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/representation.h"
#include "utils/zobrist_hasher.h"

class IPopulationCache {
public:
    virtual ~IPopulationCache() = default;
    
    // Get or calculate fitness for an individual
    virtual double getOrCalculateFitness(void* individual, 
                                       const DNAInstance& instance,
                                       std::shared_ptr<IFitness> fitness,
                                       std::shared_ptr<IRepresentation> representation) = 0;
    
    // Store a new population with their fitness values
    virtual void updatePopulation(const std::vector<void*>& population, 
                                const DNAInstance& instance,
                                std::shared_ptr<IFitness> fitness,
                                std::shared_ptr<IRepresentation> representation) = 0;
    
    // Clear the cache
    virtual void clear() = 0;
    
    // Get the size of the cache
    virtual size_t size() const = 0;
};

class CachedPopulation : public IPopulationCache {
private:
    struct CacheEntry {
        double fitness;
        uint64_t hash;
        std::chrono::steady_clock::time_point lastAccess;
    };
    
    std::unordered_map<uint64_t, CacheEntry> cache;
    const size_t maxCacheSize = 10000; // Configurable
    
    uint64_t computeHash(void* individual) {
        auto* perm = static_cast<std::vector<int>*>(individual);
        return ZobristHasher::getInstance().hashPermutation(*perm);
    }
    
    void evictOldEntries() {
        if (cache.size() <= maxCacheSize) return;
        
        std::vector<std::pair<uint64_t, std::chrono::steady_clock::time_point>> entries;
        for (const auto& entry : cache) {
            entries.emplace_back(entry.first, entry.second.lastAccess);
        }
        
        // Sort by time, oldest first
        std::sort(entries.begin(), entries.end(),
                 [](const auto& a, const auto& b) {
                     return a.second < b.second;
                 });
        
        // Remove oldest 20% of entries
        size_t toRemove = cache.size() / 5;
        for (size_t i = 0; i < toRemove; ++i) {
            cache.erase(entries[i].first);
        }
    }

public:
    double getOrCalculateFitness(void* individual, 
                                const DNAInstance& instance,
                                std::shared_ptr<IFitness> fitness,
                                std::shared_ptr<IRepresentation> representation) override {
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
    
    void updatePopulation(const std::vector<void*>& population,
                         const DNAInstance& instance,
                         std::shared_ptr<IFitness> fitness,
                         std::shared_ptr<IRepresentation> representation) override {
        for (auto* individual : population) {
            getOrCalculateFitness(individual, instance, fitness, representation);
        }
    }
    
    void clear() override {
        cache.clear();
    }
    
    size_t size() const override {
        return cache.size();
    }
};

#endif //POPULATION_CACHE_H
