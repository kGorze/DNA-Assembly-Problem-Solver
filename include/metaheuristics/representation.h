#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "../dna/dna_instance.h"
#include "../interfaces/i_representation.h"
#include "individual.h"

class PermutationRepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<Individual>> initializePopulation(size_t populationSize, const DNAInstance& instance) override;
    bool initializeIndividual(Individual& individual, const DNAInstance& instance);
    bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override;
    std::string toString(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override;
    std::string toDNA(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override;
    void clearValidationCache() const { m_validationCache.clear(); }
                        
private:
    bool validateGenes(const std::vector<int>& genes, const DNAInstance& instance) const;
    size_t countMismatches(const std::string& str1, const std::string& str2, size_t overlapLen) const;
    
    // Cache for validation results to avoid redundant validations
    mutable std::unordered_map<const Individual*, bool> m_validationCache;
}; 