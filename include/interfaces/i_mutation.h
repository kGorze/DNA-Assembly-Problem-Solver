#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "i_representation.h"

class IMutation {
public:
    virtual ~IMutation() = default;
    
    // Mutate a solution
    virtual void mutate(
        std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
}; 