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
    bool m_diversityTrackingEnabled = false;
    double m_diversityThreshold = 0.2;

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

    void enableDiversityTracking(bool enabled) override {
        m_diversityTrackingEnabled = enabled;
    }

    void setDiversityThreshold(double threshold) override {
        m_diversityThreshold = threshold;
    }

    bool isDiversityTrackingEnabled() const override {
        return m_diversityTrackingEnabled;
    }

    double getDiversityThreshold() const override {
        return m_diversityThreshold;
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
    bool m_diversityTrackingEnabled = false;
    double m_diversityThreshold = 0.2;

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

    void enableDiversityTracking(bool enabled) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_diversityTrackingEnabled = enabled;
    }

    void setDiversityThreshold(double threshold) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_diversityThreshold = threshold;
    }

    bool isDiversityTrackingEnabled() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_diversityTrackingEnabled;
    }

    double getDiversityThreshold() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_diversityThreshold;
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
        if (spectrum.empty() || genes.empty()) return 0.0;
        
        // Count how many spectrum elements are used
        std::vector<bool> used(spectrum.size(), false);
        int coveredCount = 0;
        
        for (int gene : genes) {
            if (gene >= 0 && static_cast<size_t>(gene) < spectrum.size() && !used[gene]) {
                used[gene] = true;
                coveredCount++;
            }
        }
        
        return static_cast<double>(coveredCount) / spectrum.size();
    }
    
    double calculateConnectivity(const std::vector<int>& genes, const DNAInstance& instance) {
        if (genes.size() < 2) return 0.0;
        
        int validConnections = 0;
        int totalConnections = genes.size() - 1;
        
        for (size_t i = 0; i < genes.size() - 1; ++i) {
            // Check if indices are valid
            if (genes[i] >= 0 && static_cast<size_t>(genes[i]) < instance.getSpectrum().size() &&
                genes[i+1] >= 0 && static_cast<size_t>(genes[i+1]) < instance.getSpectrum().size()) {
                // More lenient connectivity check
                if (std::abs(genes[i] - genes[i + 1]) <= instance.getDeltaK() + 3) {
                    validConnections++;
                }
            }
        }
        
        return totalConnections > 0 ? static_cast<double>(validConnections) / totalConnections : 0.0;
    }
}; 