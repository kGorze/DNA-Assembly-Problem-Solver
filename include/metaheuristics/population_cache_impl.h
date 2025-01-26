#pragma once

#include "interfaces/i_population_cache.h"
#include "metaheuristics/individual.h"
#include "dna/dna_instance.h"
#include <unordered_map>
#include <memory>

class SimplePopulationCache : public IPopulationCache {
private:
    std::unordered_map<std::shared_ptr<Individual>, double> m_cache;

public:
    void updatePopulation(const std::vector<std::shared_ptr<Individual>>& population) override {
        clear();
        for (const auto& individual : population) {
            add(individual);
        }
    }

    void clear() override {
        m_cache.clear();
    }

    void reserve(size_t size) override {
        m_cache.reserve(size);
    }

    void add(const std::shared_ptr<Individual>& individual) override {
        if (!contains(individual)) {
            m_cache[individual] = calculateFitness(individual);
        }
    }

    bool contains(const std::shared_ptr<Individual>& individual) const override {
        return m_cache.find(individual) != m_cache.end();
    }

    double getOrCalculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) override {
        auto it = m_cache.find(individual);
        if (it != m_cache.end()) {
            return it->second;
        }
        
        double fitness = calculateFitness(individual);
        m_cache[individual] = fitness;
        return fitness;
    }

private:
    double calculateFitness(const std::shared_ptr<Individual>& individual) {
        // Simple fitness calculation based on the genes
        const auto& genes = individual->getGenes();
        if (genes.empty()) {
            return 0.0;
        }
        
        // Example fitness: sum of genes divided by size
        double sum = 0.0;
        for (int gene : genes) {
            sum += gene;
        }
        return sum / genes.size();
    }
}; 