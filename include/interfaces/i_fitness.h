#pragma once

#include <memory>
#include "../dna/dna_instance.h"
#include "../metaheuristics/individual.h"

class IFitness {
public:
    virtual ~IFitness() = default;
    
    virtual double evaluate(const std::shared_ptr<Individual>& individual,
                          const DNAInstance& instance) = 0;
}; 