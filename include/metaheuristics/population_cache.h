//
// Created by konrad_guest on 09/01/2025.
// SMART

#ifndef POPULATION_CACHE_H
#define POPULATION_CACHE_H

#include <vector>
#include <unordered_map>
#include <chrono>
#include <memory>

#include "generator/dna_generator.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/representation.h"
#include "utils/zobrist_hasher.h"

class IPopulationCache {
public:
    virtual ~IPopulationCache() = default;
    
    virtual double getOrCalculateFitness(std::shared_ptr<std::vector<int>> individual, 
                                         const DNAInstance& instance,
                                         std::shared_ptr<IFitness> fitness,
                                         std::shared_ptr<IRepresentation> representation) = 0;
    
    virtual void updatePopulation(const std::vector<std::shared_ptr<std::vector<int>>>& population, 
                                  const DNAInstance& instance,
                                  std::shared_ptr<IFitness> fitness,
                                  std::shared_ptr<IRepresentation> representation) = 0;
    
    virtual void clear() = 0;
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
    const size_t maxCacheSize = 10000;
    
    uint64_t computeHash(std::shared_ptr<std::vector<int>> individual);
    void evictOldEntries();

public:
    double getOrCalculateFitness(std::shared_ptr<std::vector<int>> individual, 
                                 const DNAInstance& instance,
                                 std::shared_ptr<IFitness> fitness,
                                 std::shared_ptr<IRepresentation> representation) override;
    
    void updatePopulation(const std::vector<std::shared_ptr<std::vector<int>>>& population,
                          const DNAInstance& instance,
                          std::shared_ptr<IFitness> fitness,
                          std::shared_ptr<IRepresentation> representation) override;
    
    void clear() override;
    size_t size() const override;
};

#endif //POPULATION_CACHE_H
