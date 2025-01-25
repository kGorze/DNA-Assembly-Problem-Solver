#pragma once

#include <memory>
#include <vector>
#include "../dna/dna_instance.h"
#include "i_representation.h"
#include "../metaheuristics/individual.h"

class IFitness {
public:
    virtual ~IFitness() = default;
    
    virtual double calculateFitness(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) const = 0;
}; 