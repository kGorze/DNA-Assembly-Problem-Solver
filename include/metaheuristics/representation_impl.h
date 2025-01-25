#pragma once

#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <random>
#include <numeric>

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

    std::vector<int> generateRandomSolution(const DNAInstance& instance) const {
        std::vector<int> solution(instance.getSpectrum().size());
        std::iota(solution.begin(), solution.end(), 0); // Fill with 0,1,2,...
        
        // Ensure startIndex is at position 0
        auto startIt = std::find(solution.begin(), solution.end(), instance.getStartIndex());
        if (startIt != solution.begin()) {
            std::iter_swap(solution.begin(), startIt);
        }
        
        // Shuffle rest of the sequence
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(solution.begin() + 1, solution.end(), gen);
        
        return solution;
    }
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