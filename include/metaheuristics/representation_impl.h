#pragma once

#include "metaheuristics/i_representation.h"
#include "dna/dna_instance.h"
#include "utils/random.h"
#include <vector>
#include <memory>
#include <string>
#include <algorithm>
#include <random>
#include <numeric>
#include <mutex>
#include <sstream>

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

    bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return false;
        
        const auto& genes = individual->getGenes();
        if (genes.size() != static_cast<size_t>(instance.getN())) return false;
        
        for (const auto& gene : genes) {
            if (gene != 0 && gene != 1) return false;
        }
        
        return true;
    }

    std::string toString(const std::shared_ptr<Individual>& individual) const override {
        if (!individual) return "";
        
        std::stringstream ss;
        for (const auto& gene : individual->getGenes()) {
            ss << gene;
        }
        return ss.str();
    }

    std::string toDNA(const std::shared_ptr<Individual>& individual) const override {
        if (!individual) return "";
        
        std::stringstream ss;
        for (const auto& gene : individual->getGenes()) {
            ss << (gene == 1 ? 'A' : 'T');
        }
        return ss.str();
    }

    std::shared_ptr<Individual> generateRandomSolution() override {
        std::vector<int> genes;
        genes.reserve(m_instance.getN());
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 1);
        
        for (size_t i = 0; i < static_cast<size_t>(m_instance.getN()); ++i) {
            genes.push_back(dis(gen));
        }
        
        auto individual = std::make_shared<Individual>();
        individual->setGenes(std::move(genes));
        return individual;
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
    explicit GraphPathRepresentation(const DNAInstance& instance) : m_instance(instance) {}

    std::vector<std::shared_ptr<Individual>> initializePopulation(size_t populationSize) override {
        std::vector<std::shared_ptr<Individual>> population;
        population.reserve(populationSize);
        
        for (size_t i = 0; i < populationSize; ++i) {
            population.push_back(generateRandomSolution());
        }
        
        return population;
    }

    bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return false;
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return false;
        
        // Check if all indices are within bounds
        for (const auto& gene : genes) {
            if (gene < 0 || gene >= instance.getN()) return false;
        }
        
        // Check for duplicates
        std::vector<int> sorted = genes;
        std::sort(sorted.begin(), sorted.end());
        return std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end();
    }

    std::string toString(const std::shared_ptr<Individual>& individual) const override {
        if (!individual) return "";
        
        std::stringstream ss;
        const auto& genes = individual->getGenes();
        for (size_t i = 0; i < genes.size(); ++i) {
            if (i > 0) ss << " -> ";
            ss << genes[i];
        }
        return ss.str();
    }

    std::string toDNA(const std::shared_ptr<Individual>& individual) const override {
        if (!individual) return "";
        
        std::string dna;
        dna.reserve(m_instance.getN());
        
        std::vector<bool> visited(m_instance.getN(), false);
        const auto& genes = individual->getGenes();
        
        for (const auto& gene : genes) {
            visited[gene] = true;
            dna.push_back('A');
        }
        
        for (size_t i = 0; i < visited.size(); ++i) {
            if (!visited[i]) {
                dna.push_back('T');
            }
        }
        
        return dna;
    }

    std::shared_ptr<Individual> generateRandomSolution() override {
        std::vector<int> genes;
        genes.reserve(m_instance.getN());
        
        for (int i = 0; i < m_instance.getN(); ++i) {
            genes.push_back(i);
        }
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(genes.begin(), genes.end(), gen);
        
        auto individual = std::make_shared<Individual>();
        individual->setGenes(std::move(genes));
        return individual;
    }

private:
    const DNAInstance& m_instance;
}; 