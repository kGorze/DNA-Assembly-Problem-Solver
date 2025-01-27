#pragma once

#include <memory>
#include <vector>
#include "../dna/dna_instance.h"
#include "i_representation.h"
#include "../metaheuristics/individual.h"

class IMutation {
public:
    virtual ~IMutation() = default;
    
    // Mutate a solution
    virtual void mutate(
        std::shared_ptr<Individual>& individual,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) = 0;

    // Add mutation rate control
    virtual void setMutationRate(double rate) = 0;
    virtual double getMutationRate() const = 0;
}; 