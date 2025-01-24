#pragma once

#include "../interfaces/i_replacement.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>

class PartialReplacement : public IReplacement {
public:
    explicit PartialReplacement(double replacementRatio) : m_replacementRatio(replacementRatio) {}

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
};

class ElitistReplacement : public IReplacement {
public:
    std::vector<std::shared_ptr<std::vector<int>>> replace(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
        const std::vector<double>& populationFitness,
        const std::vector<double>& offspringFitness,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
}; 