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
    
    virtual const std::vector<std::shared_ptr<Individual>>& getCurrentPopulation() const = 0;
    
    virtual void clear() = 0;
    
    virtual void reserve(size_t size) = 0;
    
    virtual void add(const std::shared_ptr<Individual>& individual) = 0;
    
    virtual bool contains(const std::shared_ptr<Individual>& individual) const = 0;

    // Diversity tracking methods
    virtual void enableDiversityTracking(bool enabled) = 0;
    virtual void setDiversityThreshold(double threshold) = 0;
    virtual bool isDiversityTrackingEnabled() const = 0;
    virtual double getDiversityThreshold() const = 0;
}; 