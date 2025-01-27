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
        , m_stagnationCounter(0)
        , m_lastBestFitness(0.0)
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

            // Sort populations by fitness
            auto sortedPopulation = parents;
            auto sortedOffspring = offspring;
            
            std::sort(sortedPopulation.begin(), sortedPopulation.end(),
                [](const auto& a, const auto& b) {
                    return a && b && a->getFitness() > b->getFitness();
                });
                
            std::sort(sortedOffspring.begin(), sortedOffspring.end(),
                [](const auto& a, const auto& b) {
                    return a && b && a->getFitness() > b->getFitness();
                });

            // Check for stagnation
            if (!sortedPopulation.empty() && sortedPopulation[0]) {
                double currentBestFitness = sortedPopulation[0]->getFitness();
                if (std::abs(currentBestFitness - m_lastBestFitness) < 1e-6) {
                    m_stagnationCounter++;
                } else {
                    m_stagnationCounter = 0;
                }
                m_lastBestFitness = currentBestFitness;
            }

            // If stagnating, increase diversity by accepting more offspring
            if (m_stagnationCounter > 5) {
                numToReplace = static_cast<size_t>(parents.size() * 0.8); // Replace 80% of population
                LOG_INFO("Stagnation detected - increasing replacement ratio to 0.8");
                m_stagnationCounter = 0;
            }

            // Always keep the best individual from parents
            if (!sortedPopulation.empty() && sortedPopulation[0]) {
                newPopulation.push_back(std::make_shared<Individual>(*sortedPopulation[0]));
            }

            // Add offspring, including some that might be partially valid
            size_t offspringAdded = 0;
            for (const auto& individual : sortedOffspring) {
                if (offspringAdded >= numToReplace) break;
                if (individual) {
                    newPopulation.push_back(std::make_shared<Individual>(*individual));
                    offspringAdded++;
                }
            }

            // Fill remaining slots with diverse individuals from parents
            while (newPopulation.size() < parents.size() && !sortedPopulation.empty()) {
                // Take every third individual to maintain diversity
                for (size_t i = 1; i < sortedPopulation.size() && newPopulation.size() < parents.size(); i += 3) {
                    if (sortedPopulation[i]) {
                        newPopulation.push_back(std::make_shared<Individual>(*sortedPopulation[i]));
                    }
                }
                // If we still need more individuals, take the remaining best ones
                if (newPopulation.size() < parents.size()) {
                    for (size_t i = 1; i < sortedPopulation.size() && newPopulation.size() < parents.size(); i++) {
                        if (sortedPopulation[i] && std::find_if(newPopulation.begin(), newPopulation.end(),
                            [&](const auto& ind) {
                                return ind && ind->getGenes() == sortedPopulation[i]->getGenes();
                            }) == newPopulation.end()) {
                            newPopulation.push_back(std::make_shared<Individual>(*sortedPopulation[i]));
                        }
                    }
                }
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
    int m_stagnationCounter;
    double m_lastBestFitness;
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