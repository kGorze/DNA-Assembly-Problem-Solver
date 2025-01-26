#pragma once

#include "../interfaces/i_replacement.h"
#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <mutex>

class PartialReplacement : public IReplacement {
public:
    explicit PartialReplacement(double replacementRatio, std::shared_ptr<IPopulationCache> cache = nullptr) 
        : m_replacementRatio(replacementRatio)
        , m_fitnessCache(cache)
    {
        if (m_replacementRatio < 0.0 || m_replacementRatio > 1.0) {
            LOG_WARNING("Invalid replacement ratio " + std::to_string(m_replacementRatio) + " - clamping to valid range");
            m_replacementRatio = std::clamp(m_replacementRatio, 0.0, 1.0);
        }
    }
    
    std::vector<std::shared_ptr<Individual>> replace(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const std::vector<std::shared_ptr<Individual>>& offspring,
        [[maybe_unused]] const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        std::lock_guard<std::mutex> lock(m_mutex);

        if (parents.empty() || offspring.empty()) {
            LOG_WARNING("Empty population or offspring in replacement");
            return parents;
        }

        if (!representation) {
            LOG_ERROR("Null representation provided to replacement operator");
            return parents;
        }

        try {
            // Calculate number of individuals to replace
            size_t numToReplace = static_cast<size_t>(m_replacementRatio * parents.size());
            numToReplace = std::min(numToReplace, offspring.size());

            // Create a new population with the best individuals
            std::vector<std::shared_ptr<Individual>> newPopulation;
            newPopulation.reserve(parents.size());

            // Keep the best from the original population
            size_t keepFromPopulation = parents.size() - numToReplace;
            auto sortedPopulation = parents;
            std::sort(sortedPopulation.begin(), sortedPopulation.end(),
                [](const auto& a, const auto& b) {
                    return a && b && a->getFitness() > b->getFitness();
                });

            for (size_t i = 0; i < keepFromPopulation && i < sortedPopulation.size(); ++i) {
                if (sortedPopulation[i] && sortedPopulation[i]->isValid()) {
                    newPopulation.push_back(std::make_shared<Individual>(*sortedPopulation[i]));
                }
            }

            // Add the best offspring
            auto sortedOffspring = offspring;
            std::sort(sortedOffspring.begin(), sortedOffspring.end(),
                [](const auto& a, const auto& b) {
                    return a && b && a->getFitness() > b->getFitness();
                });

            for (size_t i = 0; i < numToReplace && i < sortedOffspring.size(); ++i) {
                if (sortedOffspring[i] && sortedOffspring[i]->isValid()) {
                    newPopulation.push_back(std::make_shared<Individual>(*sortedOffspring[i]));
                }
            }

            // If we don't have enough individuals, fill with the remaining best from population
            while (newPopulation.size() < parents.size() && keepFromPopulation < sortedPopulation.size()) {
                if (sortedPopulation[keepFromPopulation] && sortedPopulation[keepFromPopulation]->isValid()) {
                    newPopulation.push_back(std::make_shared<Individual>(*sortedPopulation[keepFromPopulation]));
                }
                ++keepFromPopulation;
            }

            LOG_INFO("Partial replacement complete. Population size: " + std::to_string(newPopulation.size()));
            return newPopulation;
        } catch (const std::exception& e) {
            LOG_ERROR("Error during replacement: " + std::string(e.what()));
            return parents;
        }
    }

private:
    double m_replacementRatio;
    std::shared_ptr<IPopulationCache> m_fitnessCache;
    mutable std::mutex m_mutex;
};

class ElitistReplacement : public IReplacement {
public:
    std::vector<std::shared_ptr<Individual>> replace(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const std::vector<std::shared_ptr<Individual>>& offspring,
        [[maybe_unused]] const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        std::lock_guard<std::mutex> lock(m_mutex);

        if (parents.empty() || offspring.empty()) {
            LOG_WARNING("Empty population or offspring in elitist replacement");
            return parents;
        }

        if (!representation) {
            LOG_ERROR("Null representation provided to elitist replacement operator");
            return parents;
        }

        try {
            // Find the best individual from the current population
            auto bestFromPopulation = std::max_element(parents.begin(), parents.end(),
                [](const auto& a, const auto& b) {
                    return !a || !b || a->getFitness() < b->getFitness();
                });

            // Create new population starting with the best individual
            std::vector<std::shared_ptr<Individual>> newPopulation;
            newPopulation.reserve(parents.size());

            if (bestFromPopulation != parents.end() && *bestFromPopulation && (*bestFromPopulation)->isValid()) {
                newPopulation.push_back(std::make_shared<Individual>(**bestFromPopulation));
            }

            // Sort offspring by fitness
            auto sortedOffspring = offspring;
            std::sort(sortedOffspring.begin(), sortedOffspring.end(),
                [](const auto& a, const auto& b) {
                    return a && b && a->getFitness() > b->getFitness();
                });

            // Add best offspring
            for (const auto& individual : sortedOffspring) {
                if (newPopulation.size() >= parents.size()) break;
                if (individual && individual->isValid()) {
                    newPopulation.push_back(std::make_shared<Individual>(*individual));
                }
            }

            // If we still need more individuals, take them from original population
            auto sortedPopulation = parents;
            std::sort(sortedPopulation.begin(), sortedPopulation.end(),
                [](const auto& a, const auto& b) {
                    return a && b && a->getFitness() > b->getFitness();
                });

            for (const auto& individual : sortedPopulation) {
                if (newPopulation.size() >= parents.size()) break;
                if (individual && individual->isValid()) {
                    newPopulation.push_back(std::make_shared<Individual>(*individual));
                }
            }

            LOG_INFO("Elitist replacement complete. Population size: " + std::to_string(newPopulation.size()));
            return newPopulation;
        } catch (const std::exception& e) {
            LOG_ERROR("Error during elitist replacement: " + std::string(e.what()));
            return parents;
        }
    }

private:
    mutable std::mutex m_mutex;
}; 