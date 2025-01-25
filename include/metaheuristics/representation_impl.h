#pragma once

#include "../interfaces/i_representation.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <random>
#include <numeric>
#include <mutex>

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
    std::vector<std::shared_ptr<Individual>> initializePopulation(
        int populationSize,
        const DNAInstance& instance) override;

    bool initializeIndividual(
        Individual& individual,
        const DNAInstance& instance) override;

    bool isValid(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const override;

    std::string toString(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const override;

    std::vector<char> toDNA(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const override;

private:
    mutable std::mutex m_mutex;
    mutable std::mt19937 m_rng{std::random_device{}()};

    std::vector<int> generateRandomSolution(const DNAInstance& instance) const {
        std::vector<int> solution(instance.getSpectrum().size());
        std::iota(solution.begin(), solution.end(), 0); // Fill with 0,1,2,...
        
        // Ensure startIndex is at position 0
        auto startIt = std::find(solution.begin(), solution.end(), instance.getStartIndex());
        if (startIt != solution.begin()) {
            std::iter_swap(solution.begin(), startIt);
        }
        
        // Shuffle rest of the sequence
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::shuffle(solution.begin() + 1, solution.end(), m_rng);
        }
        
        return solution;
    }

    bool validateSolution(const std::vector<int>& genes, const DNAInstance& instance) const {
        if (genes.empty()) {
            return false;
        }

        // Check if all indices are present exactly once
        std::vector<bool> used(instance.getSpectrum().size(), false);
        for (int gene : genes) {
            if (gene < 0 || gene >= (int)instance.getSpectrum().size()) {
                return false;
            }
            if (used[gene]) {
                return false;
            }
            used[gene] = true;
        }

        // Check if all indices were used
        return std::all_of(used.begin(), used.end(), [](bool v) { return v; });
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