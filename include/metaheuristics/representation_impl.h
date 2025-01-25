#pragma once

#include "../interfaces/i_representation.h"
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
#include <stdexcept>

class DirectDNARepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<Individual>> initializePopulation(int populationSize, const DNAInstance& instance) override {
        if (populationSize <= 0) {
            throw std::invalid_argument("Population size must be positive");
        }

        std::vector<std::shared_ptr<Individual>> population;
        population.reserve(populationSize);

        for (int i = 0; i < populationSize; ++i) {
            auto individual = std::make_shared<Individual>();
            std::vector<int> genes;
            // Generate random genes based on instance parameters
            genes.resize(instance.getN());
            std::generate(genes.begin(), genes.end(), []() {
                return Random::instance().getRandomInt(0, 3);  // 0-3 for A,C,G,T
            });
            individual->setGenes(genes);
            population.push_back(individual);
        }

        return population;
    }

    bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return false;
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return false;

        // Check if all genes are valid nucleotides (0-3)
        return std::all_of(genes.begin(), genes.end(), [](int gene) {
            return gene >= 0 && gene <= 3;
        });
    }

    std::string toString(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return "";
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return "";

        std::string result;
        result.reserve(genes.size());
        
        // Convert numeric genes to nucleotides
        for (int gene : genes) {
            switch (gene) {
                case 0: result += 'A'; break;
                case 1: result += 'C'; break;
                case 2: result += 'G'; break;
                case 3: result += 'T'; break;
                default: result += 'N'; break;
            }
        }
        
        return result;
    }

    std::string toDNA(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return "";
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return "";

        std::string dna;
        dna.reserve(genes.size());
        
        // Convert numeric genes to nucleotides
        for (int gene : genes) {
            switch (gene) {
                case 0: dna += 'A'; break;
                case 1: dna += 'C'; break;
                case 2: dna += 'G'; break;
                case 3: dna += 'T'; break;
                default: dna += 'N'; break;
            }
        }
        
        return dna;
    }
};

class PermutationRepresentation : public IRepresentation {
public:
    std::vector<std::shared_ptr<Individual>> initializePopulation(
        int populationSize,
        const DNAInstance& instance) override;

    bool initializeIndividual(
        Individual& individual,
        const DNAInstance& instance) {
        // Implementation here
        return true;
    }

    bool isValid(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const override;

    std::string toString(
        const std::shared_ptr<Individual>& solution,
        const DNAInstance& instance) const override;

    std::string toDNA(
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
    std::vector<std::shared_ptr<Individual>> initializePopulation(int populationSize, const DNAInstance& instance) override {
        if (populationSize <= 0) {
            throw std::invalid_argument("Population size must be positive");
        }

        std::vector<std::shared_ptr<Individual>> population;
        population.reserve(populationSize);

        for (int i = 0; i < populationSize; ++i) {
            auto individual = std::make_shared<Individual>();
            std::vector<int> genes;
            // Generate random permutation of indices
            genes.resize(instance.getSpectrum().size());
            std::iota(genes.begin(), genes.end(), 0);
            std::shuffle(genes.begin(), genes.end(), Random::instance().getGenerator());
            individual->setGenes(genes);
            population.push_back(individual);
        }

        return population;
    }

    bool isValid(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return false;
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return false;

        // Check if genes form a valid permutation
        std::vector<bool> used(genes.size(), false);
        for (int gene : genes) {
            if (gene < 0 || gene >= static_cast<int>(genes.size()) || used[gene]) {
                return false;
            }
            used[gene] = true;
        }
        return true;
    }

    std::string toString(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return "";
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return "";

        std::string result;
        result.reserve(genes.size() * 3);  // Estimate 3 chars per number
        
        for (size_t i = 0; i < genes.size(); ++i) {
            if (i > 0) result += " ";
            result += std::to_string(genes[i]);
        }
        
        return result;
    }

    std::string toDNA(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        if (!individual) return "";
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return "";

        const auto& spectrum = instance.getSpectrum();
        if (spectrum.empty()) return "";

        // Convert path through spectrum to DNA sequence
        std::string dna;
        dna.reserve(instance.getN());  // Reserve estimated size
        
        // Add first oligonucleotide completely
        const auto& first = spectrum[genes[0]];
        dna += first;
        
        // For subsequent oligonucleotides, add only non-overlapping part
        for (size_t i = 1; i < genes.size(); ++i) {
            const auto& oligo = spectrum[genes[i]];
            if (oligo.size() > instance.getK()) {
                dna += oligo.substr(instance.getK());
            }
        }
        
        return dna;
    }
}; 