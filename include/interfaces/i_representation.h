#pragma once

#include "metaheuristics/individual.h"
#include "dna/dna_instance.h"
#include <memory>
#include <vector>
#include <string>

class IRepresentation {
public:
    virtual ~IRepresentation() = default;

    virtual std::vector<std::shared_ptr<Individual>> initializePopulation(int populationSize, const DNAInstance& instance) = 0;
    virtual bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const = 0;
    virtual std::string toString(const std::shared_ptr<Individual>& individual) const = 0;
    virtual std::string toDNA(const std::shared_ptr<Individual>& individual) const = 0;
}; 