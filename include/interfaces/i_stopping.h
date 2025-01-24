#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"

class IStopping {
public:
    virtual ~IStopping() = default;
    
    // Check if algorithm should stop
    virtual bool stop(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        int currentGeneration,
        double bestFitness
    ) const = 0;
    
    // Reset stopping criteria
    virtual void reset() = 0;
}; 