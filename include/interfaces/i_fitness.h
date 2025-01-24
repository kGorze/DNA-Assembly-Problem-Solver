#pragma once

#include <memory>
#include <vector>
#include "../dna/dna_instance.h"
#include "i_representation.h"

class IFitness {
public:
    virtual ~IFitness() = default;
    
    // Calculate fitness for a solution
    virtual double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) const = 0;
}; 