#pragma once

#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>
#include <string>

class DirectDNARepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<std::vector<int>>> initializePopulation(
        int popSize, 
        const DNAInstance& instance) override;

    bool isValid(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

    std::string toString(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

    std::vector<char> toDNA(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

private:
    static constexpr int fallbackN = 100;  // Default length if instance.getN() <= 0
};

class PermutationRepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<std::vector<int>>> initializePopulation(
        int populationSize,
        const DNAInstance& instance) override;

    bool isValid(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

    std::string toString(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

    std::vector<char> toDNA(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;
};

class GraphPathRepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<std::vector<int>>> initializePopulation(
        int popSize,
        const DNAInstance& instance) override;

    bool isValid(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

    std::string toString(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;

    std::vector<char> toDNA(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance) const override;
}; 