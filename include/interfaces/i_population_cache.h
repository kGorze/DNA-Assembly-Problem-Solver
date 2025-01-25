#pragma once

#include <memory>
#include <vector>
#include "../dna/dna_instance.h"
#include "../metaheuristics/individual.h"

class IPopulationCache {
public:
    virtual ~IPopulationCache() = default;
    
    virtual double getOrCalculateFitness(const std::shared_ptr<Individual>& individual,
                                       const DNAInstance& instance) = 0;
    
    virtual void updatePopulation(const std::vector<std::shared_ptr<Individual>>& population) = 0;
    
    virtual void clear() = 0;
}; 