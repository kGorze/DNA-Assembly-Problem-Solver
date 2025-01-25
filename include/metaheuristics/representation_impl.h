#pragma once

#include "interfaces/i_representation.h"
#include "dna/dna_instance.h"
#include "utils/random.h"
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <random>
#include <numeric>
#include <mutex>

class DirectDNARepresentation : public IRepresentation {
public:
    explicit DirectDNARepresentation(const DNAInstance& instance) : m_instance(instance) {}

    std::vector<std::shared_ptr<Individual>> initializePopulation(size_t populationSize) override {
        std::vector<std::shared_ptr<Individual>> population;
        population.reserve(populationSize);
        
        for (size_t i = 0; i < populationSize; ++i) {
            population.push_back(generateRandomSolution());
        }
        
        return population;
    }

    bool isValid(const std::shared_ptr<Individual>& individual) const override {
        if (!individual || individual->empty()) return false;
        
        const auto& genes = individual->getGenes();
        if (genes.size() != m_instance.getN()) return false;
        
        // Check if all values are valid nucleotides (0-3)
        return std::all_of(genes.begin(), genes.end(), 
            [](int gene) { return gene >= 0 && gene <= 3; });
    }

    std::string toString(const std::shared_ptr<Individual>& individual) const override {
        if (!individual) return "";
        
        const auto& genes = individual->getGenes();
        std::string result;
        result.reserve(genes.size());
        
        for (int gene : genes) {
            switch (gene) {
                case 0: result += 'A'; break;
                case 1: result += 'C'; break;
                case 2: result += 'G'; break;
                case 3: result += 'T'; break;
                default: result += 'X'; break;
            }
        }
        
        return result;
    }

    std::shared_ptr<Individual> generateRandomSolution() override {
        auto individual = std::make_shared<Individual>();
        std::vector<int> genes;
        genes.reserve(m_instance.getN());
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 3);
        
        for (size_t i = 0; i < m_instance.getN(); ++i) {
            genes.push_back(dis(gen));
        }
        
        individual->setGenes(std::move(genes));
        individual->setValid(true);
        return individual;
    }

    std::string toDNA(const std::shared_ptr<Individual>& individual) const override {
        return toString(individual);
    }

private:
    const DNAInstance& m_instance;
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
    std::vector<std::shared_ptr<Individual>> initializePopulation(
        int populationSize,
        const DNAInstance& instance) override {
        std::vector<std::shared_ptr<Individual>> population;
        population.reserve(populationSize);
        
        for (int i = 0; i < populationSize; ++i) {
            auto individual = std::make_shared<Individual>();
            individual->setGenes(generateRandomSolution(instance));
            population.push_back(individual);
        }
        
        return population;
    }

    bool isValid(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance) const override {
        if (!individual || individual->empty()) return false;
        
        const auto& genes = individual->getGenes();
        const auto& spectrum = instance.getSpectrum();
        
        // Check if all indices are within bounds
        return std::all_of(genes.begin(), genes.end(),
            [&spectrum](int gene) { return gene >= 0 && gene < static_cast<int>(spectrum.size()); });
    }

    std::string toString(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance) const override {
        if (!individual || individual->empty()) return "";
        
        const auto& genes = individual->getGenes();
        const auto& spectrum = instance.getSpectrum();
        
        std::string result;
        for (int gene : genes) {
            if (gene >= 0 && gene < static_cast<int>(spectrum.size())) {
                result += spectrum[gene];
            }
        }
        
        return result;
    }

    std::vector<int> generateRandomSolution(
        const DNAInstance& instance) override {
        const auto& spectrum = instance.getSpectrum();
        if (spectrum.empty()) return {};
        
        std::vector<int> solution;
        solution.reserve(spectrum.size());
        
        // Create a permutation of spectrum indices
        for (size_t i = 0; i < spectrum.size(); ++i) {
            solution.push_back(i);
        }
        
        // Shuffle the indices
        Random random;
        std::shuffle(solution.begin(), solution.end(), random.getGenerator());
        
        return solution;
    }

    std::vector<char> toDNA(
        const std::shared_ptr<Individual>& individual,
        const DNAInstance& instance) const override {
        if (!individual || individual->empty()) return {};
        
        const auto& genes = individual->getGenes();
        const auto& spectrum = instance.getSpectrum();
        
        std::vector<char> dna;
        for (int gene : genes) {
            if (gene >= 0 && gene < static_cast<int>(spectrum.size())) {
                const auto& oligo = spectrum[gene];
                dna.insert(dna.end(), oligo.begin(), oligo.end());
            }
        }
        
        return dna;
    }
}; 