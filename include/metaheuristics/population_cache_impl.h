#pragma once

#include <unordered_map>
#include <string>
#include "../dna/dna_instance.h"
#include "../interfaces/i_population_cache.h"
#include "metaheuristics/individual.h"
#include <memory>
#include <vector>
#include <mutex>

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

    double getOrCalculateFitness(const std::shared_ptr<Individual>& individual, 
                               [[maybe_unused]] const DNAInstance& instance) override {
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

class PopulationCache : public IPopulationCache {
private:
    std::vector<std::shared_ptr<Individual>> m_population;
    std::unordered_map<std::shared_ptr<Individual>, double> m_fitnessCache;
    mutable std::mutex m_mutex;

public:
    void updatePopulation(const std::vector<std::shared_ptr<Individual>>& population) override {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_population = population;
        m_fitnessCache.clear();
    }
    
    const std::vector<std::shared_ptr<Individual>>& getCurrentPopulation() const override {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_population;
    }
    
    void clear() override {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_population.clear();
        m_fitnessCache.clear();
    }

    void reserve(size_t size) override {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_population.reserve(size);
        m_fitnessCache.reserve(size);
    }

    void add(const std::shared_ptr<Individual>& individual) override {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (!contains(individual)) {
            m_population.push_back(individual);
        }
    }

    bool contains(const std::shared_ptr<Individual>& individual) const override {
        std::lock_guard<std::mutex> lock(m_mutex);
        return std::find(m_population.begin(), m_population.end(), individual) != m_population.end();
    }

    double getOrCalculateFitness(const std::shared_ptr<Individual>& individual,
                                const DNAInstance& instance) override {
        std::lock_guard<std::mutex> lock(m_mutex);
        auto it = m_fitnessCache.find(individual);
        if (it != m_fitnessCache.end()) {
            return it->second;
        }
        
        // Calculate fitness using instance parameters
        double fitness = calculateFitness(individual, instance);
        m_fitnessCache[individual] = fitness;
        return fitness;
    }

private:
    double calculateFitness(const std::shared_ptr<Individual>& individual,
                          const DNAInstance& instance) {
        if (!individual) return 0.0;
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return 0.0;
        
        // Calculate fitness based on DNA sequence properties
        double coverage = calculateCoverage(genes, instance);
        double connectivity = calculateConnectivity(genes, instance);
        
        return (coverage + connectivity) / 2.0;
    }
    
    double calculateCoverage(const std::vector<int>& genes, const DNAInstance& instance) {
        const auto& spectrum = instance.getSpectrum();
        int coveredCount = 0;
        
        for (const auto& target : spectrum) {
            bool found = false;
            for (size_t i = 0; i <= genes.size() - target.size(); ++i) {
                bool match = true;
                for (size_t j = 0; j < target.size(); ++j) {
                    if (genes[i + j] != target[j]) {
                        match = false;
                        break;
                    }
                }
                if (match) {
                    found = true;
                    break;
                }
            }
            if (found) coveredCount++;
        }
        
        return spectrum.empty() ? 0.0 : static_cast<double>(coveredCount) / spectrum.size();
    }
    
    double calculateConnectivity(const std::vector<int>& genes, const DNAInstance& instance) {
        if (genes.size() < 2) return 0.0;
        
        int validConnections = 0;
        for (size_t i = 0; i < genes.size() - 1; ++i) {
            if (std::abs(genes[i] - genes[i + 1]) <= instance.getDeltaK()) {
                validConnections++;
            }
        }
        
        return static_cast<double>(validConnections) / (genes.size() - 1);
    }
}; 