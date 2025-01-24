#pragma once

#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include <unordered_map>
#include <shared_mutex>
#include <memory>
#include <sstream>
#include <string>

class SimplePopulationCache : public IPopulationCache {
public:
    static constexpr size_t MAX_CACHE_SIZE = 10000;

    SimplePopulationCache() = default;

    double getOrCalculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<const IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) override;

    void updatePopulation(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) override;

    void clear() override;

    size_t size() const {
        std::shared_lock<std::shared_mutex> lock(m_mutex);
        return m_cache.size();
    }

private:
    std::unordered_map<std::string, double> m_cache;
    mutable std::shared_mutex m_mutex;
    
    std::string getSolutionKey(const std::vector<int>& solution);
}; 