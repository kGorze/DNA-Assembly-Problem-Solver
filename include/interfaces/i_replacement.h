#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "i_fitness.h"
#include "i_representation.h"

class IReplacement {
public:
    virtual ~IReplacement() = default;
    
    // Replace old population with offspring
    virtual std::vector<std::shared_ptr<std::vector<int>>> replace(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
}; 