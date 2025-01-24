#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"

class IRepresentation {
public:
    virtual ~IRepresentation() = default;
    
    // Initialize a population of solutions
    virtual std::vector<std::shared_ptr<std::vector<int>>> initializePopulation(
        int populationSize, 
        const DNAInstance& instance
    ) = 0;
    
    // Validate if a solution is valid for the given instance
    virtual bool isValid(
        const std::shared_ptr<std::vector<int>>& solution, 
        const DNAInstance& instance
    ) const = 0;
}; 