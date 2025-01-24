#pragma once

#include "../interfaces/i_fitness.h"
#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>

class SimpleFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
};

class BetterFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
};

class SmithWatermanFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

private:
    double smithWaterman(const std::string& seq1, const std::string& seq2);
};

class OptimizedGraphBasedFitness : public IFitness {
public:
    double calculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
}; 