#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "i_fitness.h"
#include "i_representation.h"

class ISelection {
public:
    virtual ~ISelection() = default;
    
    // Select parents from the population for breeding
    virtual std::vector<std::shared_ptr<std::vector<int>>> select(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
}; 