#pragma once

#include <memory>
#include <vector>
#include <string>
#include "../dna/dna_instance.h"
#include "../interfaces/i_representation.h"
#include "individual.h"

class PermutationRepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<Individual>> initializePopulation(int populationSize, const DNAInstance& instance) override;
    bool initializeIndividual(Individual& individual, const DNAInstance& instance);
    bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override;
    std::string toString(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override;
    std::string toDNA(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override;
                        
private:
    bool validateGenes(const std::vector<int>& genes, const DNAInstance& instance) const;
}; 