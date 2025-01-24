#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "i_fitness.h"
#include "i_representation.h"

class IPopulationCache {
public:
    virtual ~IPopulationCache() = default;
    
    // Get or calculate fitness for a solution
    virtual double getOrCalculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
    
    // Update cache with new population
    virtual void updatePopulation(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
    
    // Clear the cache
    virtual void clear() = 0;
}; 