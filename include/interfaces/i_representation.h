#pragma once

#include <memory>
#include <string>
#include <vector>
#include "../metaheuristics/individual.h"
#include "../dna/dna_instance.h"

class IRepresentation {
public:
    virtual ~IRepresentation() = default;
    
    virtual bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const = 0;
    virtual std::string toString(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const = 0;
    virtual std::string toDNA(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const = 0;
    virtual std::vector<std::shared_ptr<Individual>> initializePopulation(size_t populationSize, const DNAInstance& instance) = 0;
    virtual bool initializeIndividual(Individual& individual, const DNAInstance& instance) = 0;
}; 