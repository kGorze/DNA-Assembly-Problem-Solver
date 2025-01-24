#pragma once

#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include <unordered_map>
#include <mutex>
#include <memory>
#include <sstream>

class SimplePopulationCache : public IPopulationCache {
public:
    double getOrCalculateFitness(
        const std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitnessCalculator,
        std::shared_ptr<IRepresentation> representation
    ) override;

    void updatePopulation(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) override {
        std::lock_guard<std::mutex> lock(cacheMutex);
        fitnessCache.clear();
        for (const auto& individual : population) {
            if (individual) {
                std::stringstream ss;
                for (int gene : *individual) {
                    ss << gene << ",";
                }
                fitnessCache[ss.str()] = fitness->calculateFitness(individual, instance, representation);
            }
        }
    }

    void clear() override;

private:
    std::unordered_map<std::string, double> fitnessCache;
    std::mutex cacheMutex;
    
    std::string getSolutionKey(const std::vector<int>& solution);
}; 