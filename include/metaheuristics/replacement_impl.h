#pragma once

#include "../interfaces/i_replacement.h"
#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>

class PartialReplacement : public IReplacement {
public:
    explicit PartialReplacement(double replacementRatio, std::shared_ptr<IPopulationCache> cache = nullptr) 
        : m_replacementRatio(replacementRatio), m_fitnessCache(cache) {}
    
    virtual ~PartialReplacement() = default;
    
    std::vector<std::shared_ptr<std::vector<int>>> replace(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
        const std::vector<double>& populationFitness,
        const std::vector<double>& offspringFitness,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

private:
    double m_replacementRatio;
    std::shared_ptr<IPopulationCache> m_fitnessCache;
};

class ElitistReplacement : public IReplacement {
public:
    virtual ~ElitistReplacement() = default;
    
    std::vector<std::shared_ptr<std::vector<int>>> replace(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
        const std::vector<double>& populationFitness,
        const std::vector<double>& offspringFitness,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
}; 