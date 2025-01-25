#pragma once

#include <vector>
#include <memory>
#include <string>
#include "../dna/dna_instance.h"
#include "../metaheuristics/individual.h"

class IRepresentation {
public:
    virtual ~IRepresentation() = default;
    
    // Initialize a population of solutions
    virtual std::vector<std::shared_ptr<Individual>> initializePopulation(
        int populationSize, 
        const DNAInstance& instance
    ) = 0;
    
    // Initialize a single individual
    virtual bool initializeIndividual(
        Individual& individual,
        const DNAInstance& instance
    ) = 0;
    
    // Validate if a solution is valid for the given instance
    virtual bool isValid(
        const std::shared_ptr<Individual>& solution, 
        const DNAInstance& instance
    ) const = 0;

    // Convert solution to string representation
    virtual std::string toString(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance
    ) const = 0;

    // Convert solution to DNA sequence
    virtual std::vector<char> toDNA(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance
    ) const = 0;
}; 