#pragma once

#include <vector>
#include <memory>
#include "../dna/dna_instance.h"
#include "i_representation.h"

class IReplacement {
public:
    virtual ~IReplacement() = default;
    
    virtual std::vector<std::shared_ptr<std::vector<int>>> replace(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
        const std::vector<double>& populationFitness,
        const std::vector<double>& offspringFitness,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
}; 