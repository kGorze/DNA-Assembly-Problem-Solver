#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "../metaheuristics/individual.h"

class IStopping {
public:
    virtual ~IStopping() = default;
    
    // Check if algorithm should stop
    virtual bool shouldStop(
        int currentGeneration,
        double bestFitness
    ) const = 0;
    
    // Reset stopping criteria
    virtual void reset() = 0;
}; 