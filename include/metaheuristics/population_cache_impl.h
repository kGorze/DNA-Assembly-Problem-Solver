#pragma once

#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include <unordered_map>
#include <mutex>
#include <memory>
#include <sstream>

class SimplePopulationCache : public IPopulationCache {
public:
    double getOrCalculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitnessCalculator,
        std::shared_ptr<IRepresentation> representation
    ) override;

    void clear() override;

private:
    std::unordered_map<std::string, double> m_cache;
    mutable std::mutex m_mutex;
    
    std::string getSolutionKey(const std::vector<int>& solution);
}; 