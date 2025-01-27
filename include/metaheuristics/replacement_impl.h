#pragma once

#include "../interfaces/i_replacement.h"
#include "../interfaces/i_population_cache.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include <vector>
#include <memory>
#include <algorithm>
#include <mutex>
#include <cmath>

// Helper function to calculate distance between individuals
inline double calculateDistance(const std::shared_ptr<Individual>& ind1, 
                              const std::shared_ptr<Individual>& ind2) {
    if (!ind1 || !ind2) return std::numeric_limits<double>::max();
    
    const auto& genes1 = ind1->getGenes();
    const auto& genes2 = ind2->getGenes();
    
    if (genes1.size() != genes2.size()) return std::numeric_limits<double>::max();
    
    // Calculate Hamming distance
    int distance = 0;
    for (size_t i = 0; i < genes1.size(); ++i) {
        if (genes1[i] != genes2[i]) distance++;
    }
    
    return static_cast<double>(distance) / genes1.size();
}

// Helper function to calculate sharing factor
inline double calculateSharingFactor(double distance, double sharingRadius) {
    if (distance >= sharingRadius) return 0.0;
    return 1.0 - std::pow(distance / sharingRadius, 2);
}

class PartialReplacement : public IReplacement {
public:
    explicit PartialReplacement(double replacementRatio = 0.8, 
                              double sharingRadius = 0.2,
                              std::shared_ptr<IPopulationCache> cache = nullptr) 
        : m_replacementRatio(replacementRatio)
        , m_sharingRadius(sharingRadius)
        , m_fitnessCache(cache)
        , m_stagnationCounter(0)
        , m_lastBestFitness(0.0)
    {
        if (m_replacementRatio < 0.0 || m_replacementRatio > 1.0) {
            LOG_WARNING("Invalid replacement ratio " + std::to_string(m_replacementRatio) + 
                       " - clamping to valid range");
            m_replacementRatio = std::clamp(m_replacementRatio, 0.0, 1.0);
        }
    }
    
    std::vector<std::shared_ptr<Individual>> replace(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const std::vector<std::shared_ptr<Individual>>& offspring,
        const DNAInstance& instance,
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

            // Create a new population starting with the best parent
            std::vector<std::shared_ptr<Individual>> newPopulation;
            newPopulation.reserve(parents.size());

            // Apply fitness sharing to both populations
            auto sharedFitnessParents = calculateSharedFitness(parents);
            auto sharedFitnessOffspring = calculateSharedFitness(offspring);

            // Sort populations by shared fitness
            auto sortedPopulation = parents;
            auto sortedOffspring = offspring;
            
            std::sort(sortedPopulation.begin(), sortedPopulation.end(),
                [&](const auto& a, const auto& b) {
                    return sharedFitnessParents[a.get()] > sharedFitnessParents[b.get()];
                });
                
            std::sort(sortedOffspring.begin(), sortedOffspring.end(),
                [&](const auto& a, const auto& b) {
                    return sharedFitnessOffspring[a.get()] > sharedFitnessOffspring[b.get()];
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
                numToReplace = static_cast<size_t>(parents.size() * 0.9); // Replace 90% of population
                LOG_INFO("Stagnation detected - increasing replacement ratio to 0.9");
                m_stagnationCounter = 0;
            }

            // Always keep the best individual from parents
            if (!sortedPopulation.empty() && sortedPopulation[0]) {
                newPopulation.push_back(std::make_shared<Individual>(*sortedPopulation[0]));
            }

            // Apply crowding: for each offspring, replace the most similar parent
            std::vector<bool> parentUsed(parents.size(), false);
            parentUsed[0] = true;  // Mark best parent as used

            for (const auto& offspring : sortedOffspring) {
                if (newPopulation.size() >= parents.size()) break;
                
                // Find most similar parent that hasn't been used
                size_t mostSimilarIdx = 1;  // Start from 1 to preserve best parent
                double minDistance = std::numeric_limits<double>::max();
                
                for (size_t i = 1; i < parents.size(); ++i) {
                    if (!parentUsed[i]) {
                        double distance = calculateDistance(offspring, parents[i]);
                        if (distance < minDistance) {
                            minDistance = distance;
                            mostSimilarIdx = i;
                        }
                    }
                }
                
                // Replace the most similar parent if offspring is better
                if (!parentUsed[mostSimilarIdx] && 
                    sharedFitnessOffspring[offspring.get()] > 
                    sharedFitnessParents[parents[mostSimilarIdx].get()]) {
                    newPopulation.push_back(std::make_shared<Individual>(*offspring));
                    parentUsed[mostSimilarIdx] = true;
                } else {
                    newPopulation.push_back(std::make_shared<Individual>(*parents[mostSimilarIdx]));
                    parentUsed[mostSimilarIdx] = true;
                }
            }

            // Fill remaining slots with unused parents
            for (size_t i = 1; i < parents.size() && newPopulation.size() < parents.size(); ++i) {
                if (!parentUsed[i]) {
                    newPopulation.push_back(std::make_shared<Individual>(*parents[i]));
                    parentUsed[i] = true;
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
    std::unordered_map<Individual*, double> calculateSharedFitness(
        const std::vector<std::shared_ptr<Individual>>& population) {
        std::unordered_map<Individual*, double> sharedFitness;
        
        for (const auto& ind1 : population) {
            if (!ind1) continue;
            
            double niche_count = 0.0;
            for (const auto& ind2 : population) {
                if (!ind2) continue;
                
                double distance = calculateDistance(ind1, ind2);
                niche_count += calculateSharingFactor(distance, m_sharingRadius);
            }
            
            sharedFitness[ind1.get()] = ind1->getFitness() / std::max(1.0, niche_count);
        }
        
        return sharedFitness;
    }

    double m_replacementRatio;
    double m_sharingRadius;
    std::shared_ptr<IPopulationCache> m_fitnessCache;
    mutable std::mutex m_mutex;
    int m_stagnationCounter;
    double m_lastBestFitness;
};

class GenerationalReplacement : public IReplacement {
public:
    std::vector<std::shared_ptr<Individual>> replace(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const std::vector<std::shared_ptr<Individual>>& offspring,
        [[maybe_unused]] const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) override {
        std::lock_guard<std::mutex> lock(m_mutex);

        if (parents.empty() || offspring.empty()) {
            LOG_WARNING("Empty population or offspring in generational replacement");
            return parents;
        }

        if (!representation) {
            LOG_ERROR("Null representation provided to generational replacement operator");
            return parents;
        }

        try {
            // Find the best individual from the current population
            auto bestFromParents = std::max_element(parents.begin(), parents.end(),
                [](const auto& a, const auto& b) {
                    return !a || !b || a->getFitness() < b->getFitness();
                });

            // Create new population starting with offspring
            std::vector<std::shared_ptr<Individual>> newPopulation;
            newPopulation.reserve(parents.size());

            // Add the best parent (elitism)
            if (bestFromParents != parents.end() && *bestFromParents) {
                newPopulation.push_back(std::make_shared<Individual>(**bestFromParents));
            }

            // Add offspring
            for (const auto& individual : offspring) {
                if (newPopulation.size() >= parents.size()) break;
                if (individual) {
                    newPopulation.push_back(std::make_shared<Individual>(*individual));
                }
            }

            // If we need more individuals, fill with best parents
            if (newPopulation.size() < parents.size()) {
                auto sortedParents = parents;
                std::sort(sortedParents.begin(), sortedParents.end(),
                    [](const auto& a, const auto& b) {
                        return a && b && a->getFitness() > b->getFitness();
                    });

                for (size_t i = 1; i < sortedParents.size() && newPopulation.size() < parents.size(); ++i) {
                    if (sortedParents[i]) {
                        newPopulation.push_back(std::make_shared<Individual>(*sortedParents[i]));
                    }
                }
            }

            LOG_INFO("Generational replacement complete. Population size: " + 
                    std::to_string(newPopulation.size()));
            return newPopulation;
        } catch (const std::exception& e) {
            LOG_ERROR("Error during generational replacement: " + std::string(e.what()));
            return parents;
        }
    }

private:
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