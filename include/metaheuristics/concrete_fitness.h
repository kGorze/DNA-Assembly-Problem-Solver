#pragma once

#include "fitness_impl.h"
#include <memory>

class ConcreteOptimizedFitness : public OptimizedGraphBasedFitness {
public:
    explicit ConcreteOptimizedFitness(const GAConfig& config) 
        : OptimizedGraphBasedFitness(config) {}

    ConcreteOptimizedFitness(const GAConfig& config,
                            const std::vector<std::shared_ptr<Individual>>& population)
        : OptimizedGraphBasedFitness(config) {
        setCurrentPopulation(population);
    }

    const std::vector<std::shared_ptr<Individual>>& getCurrentPopulation() const {
        return OptimizedGraphBasedFitness::getCurrentPopulation();
    }

    double calculateEdgeQuality(const std::shared_ptr<Individual>& individual,
                              const DNAInstance& instance) const {
        return OptimizedGraphBasedFitness::calculateEdgeQuality(individual, instance);
    }

    double calculateLength(const std::shared_ptr<Individual>& individual,
                         const DNAInstance& instance) const {
        return OptimizedGraphBasedFitness::calculateLength(individual, instance);
    }

    double calculateCoverage(const std::shared_ptr<Individual>& individual,
                           const DNAInstance& instance) const {
        if (!individual) return 0.0;
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return 0.0;
        
        // Calculate coverage based on spectrum coverage
        int coveredSpectra = 0;
        const auto& spectrum = instance.getSpectrum();
        
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
            if (found) coveredSpectra++;
        }
        
        return static_cast<double>(coveredSpectra) / spectrum.size();
    }

    double calculateConnectivity(const std::shared_ptr<Individual>& individual,
                               const DNAInstance& instance) const {
        if (!individual) return 0.0;
        
        const auto& genes = individual->getGenes();
        if (genes.empty()) return 0.0;
        
        // Calculate connectivity based on overlaps between adjacent k-mers
        int validConnections = 0;
        int totalConnections = genes.size() - instance.getK() + 1;
        
        for (size_t i = 0; i < genes.size() - instance.getK(); ++i) {
            bool validConnection = true;
            for (int j = 1; j < instance.getK(); ++j) {
                if (genes[i + j] != genes[i + j + 1]) {
                    validConnection = false;
                    break;
                }
            }
            if (validConnection) validConnections++;
        }
        
        return totalConnections > 0 ? 
               static_cast<double>(validConnections) / totalConnections : 0.0;
    }
}; 