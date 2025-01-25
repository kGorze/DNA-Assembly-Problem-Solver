#pragma once

#include <memory>
#include <vector>
#include "../dna/dna_instance.h"
#include "../metaheuristics/individual.h"
#include "i_representation.h"

class ICrossover {
public:
    virtual ~ICrossover() = default;
    
    // Perform crossover between parents to create offspring
    virtual std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) = 0;
}; 