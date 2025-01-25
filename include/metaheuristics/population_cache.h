#pragma once

#include "interfaces/i_population_cache.h"
#include "interfaces/i_fitness.h"
#include "../dna/dna_instance.h"
#include "individual.h"
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <mutex>

class SimplePopulationCache : public IPopulationCache {
private:
    static constexpr size_t maxCacheSize = 10000;
    std::unordered_map<std::string, double> cache;
    std::mutex cacheMutex;
    std::shared_ptr<IFitness> m_fitness;
    
    void cleanupCache();
    
public:
    explicit SimplePopulationCache(std::shared_ptr<IFitness> fitness)
        : m_fitness(std::move(fitness)) {}
    ~SimplePopulationCache() override = default;
    
    double getOrCalculateFitness(const std::shared_ptr<Individual>& individual,
                                const DNAInstance& instance) override;
    
    void updatePopulation(const std::vector<std::shared_ptr<Individual>>& population) override;
    
    void clear() override;
}; 