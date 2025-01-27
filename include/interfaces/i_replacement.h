#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "i_representation.h"
#include "../metaheuristics/individual.h"

class IReplacement {
public:
    virtual ~IReplacement() = default;
    
    // Replace old population with offspring
    virtual std::vector<std::shared_ptr<Individual>> replace(
        const std::vector<std::shared_ptr<Individual>>& population,
        const std::vector<std::shared_ptr<Individual>>& offspring,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) = 0;

    // Add replacement ratio control
    virtual void setReplacementRatio(double ratio) = 0;
    virtual double getReplacementRatio() const = 0;
}; 