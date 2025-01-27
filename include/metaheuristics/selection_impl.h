#pragma once

#include "../interfaces/i_selection.h"
#include "../interfaces/i_population_cache.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include "utils/random.h"
#include "individual.h"
#include <vector>
#include <memory>
#include <random>
#include <mutex>
#include <algorithm>

class TournamentSelection : public ISelection {
public:
    TournamentSelection(const GAConfig& config, std::shared_ptr<IPopulationCache> cache = nullptr);

    std::vector<std::shared_ptr<Individual>> select(
        const std::vector<std::shared_ptr<Individual>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation) override;

private:
    const GAConfig& m_config;
    std::shared_ptr<IPopulationCache> m_cache;
};

class RankBasedSelection : public ISelection {
private:
    const GAConfig& m_config;
    
    // Helper function to calculate selection probabilities based on rank
    std::vector<double> calculateRankProbabilities(int populationSize) const {
        std::vector<double> probs(populationSize);
        double sum = populationSize * (populationSize + 1) / 2.0;
        for (int i = 0; i < populationSize; i++) {
            probs[i] = (populationSize - i) / sum;
        }
        return probs;
    }

public:
    explicit RankBasedSelection(const GAConfig& config) : m_config(config) {}

    std::vector<std::shared_ptr<Individual>> select(
        const std::vector<std::shared_ptr<Individual>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation) override {
        
        if (population.empty()) return {};

        // Sort population by fitness
        std::vector<std::pair<std::shared_ptr<Individual>, double>> sortedPop;
        sortedPop.reserve(population.size());
        
        for (const auto& ind : population) {
            if (!ind) continue;
            double fit = fitness->calculateFitness(ind, instance, representation);
            sortedPop.emplace_back(ind, fit);
        }
        
        std::sort(sortedPop.begin(), sortedPop.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });

        // Calculate selection probabilities based on rank
        auto probs = calculateRankProbabilities(sortedPop.size());
        
        // Select individuals using roulette wheel selection
        std::vector<std::shared_ptr<Individual>> selected;
        selected.reserve(population.size());
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        while (selected.size() < population.size()) {
            double r = dis(gen);
            double sum = 0.0;
            
            for (size_t i = 0; i < sortedPop.size(); i++) {
                sum += probs[i];
                if (r <= sum) {
                    selected.push_back(std::make_shared<Individual>(*sortedPop[i].first));
                    break;
                }
            }
        }
        
        return selected;
    }
}; 