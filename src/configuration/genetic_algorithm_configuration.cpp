//
// Created by konrad_guest on 23/01/2025.
//
#include "../../include/configuration/genetic_algorithm_configuration.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
#include "../../include/metaheuristics/concrete_fitness.h"
#include "../../include/metaheuristics/population_cache_impl.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <mutex>
#include <shared_mutex>
#include <limits>
#include "../../include/metaheuristics/selection_impl.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/replacement_impl.h"
#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
#include "../../include/interfaces/i_population_cache.h"
#include "../../include/metaheuristics/population_cache_impl.h"
#include "../../include/metaheuristics/representation.h"
#include "../../include/metaheuristics/adaptive_crossover.h"
#include "../../include/utils/logging.h"
#include <stdexcept>
#include <string>

// Add mutex declaration at the top level
static std::shared_mutex configMutex;

/*
    Konstruktor prywatny: ustawiamy domyślne wartości.
*/
// GAConfig::GAConfig() { ... }

void GAConfig::resetToDefaults()
{
    m_populationSize = 100;
    m_mutationRate = 0.2;
    m_crossoverProbability = 0.8;
    m_targetFitness = 1.0;
    m_tournamentSize = 3;  // Changed from 1 to 3 for better selection pressure
    m_replacementRatio = 0.7;
    m_selectionMethod = "rank";
    m_noImprovementGenerations = 30;
    m_timeLimitSeconds = 60;
    
    // Instance-specific parameters with defaults
    m_k = 7;
    m_deltaK = 0;
    m_lNeg = 10;
    m_lPoz = 10;
    m_repAllowed = true;
    m_probablePositive = 0;
    
    // Reset adaptive parameters
    m_adaptiveParams.useAdaptiveMutation = true;
    m_adaptiveParams.minMutationRate = 0.1;
    m_adaptiveParams.maxMutationRate = 0.4;
    m_adaptiveParams.stagnationGenerations = 5;
    m_adaptiveParams.improvementThreshold = 0.01;
    
    // Reset diversity parameters
    m_diversityParams.useFitnessSharing = true;
    m_diversityParams.useCrowding = false;
    m_diversityParams.sharingRadius = 0.2;
    m_diversityParams.sharingAlpha = 1.0;
    m_diversityParams.diversityWeight = 0.3;
}

bool GAConfig::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        LOG_ERROR("Failed to open configuration file: " + filename);
        return false;
    }

    // Start with default values
    resetToDefaults();
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            // Trim whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);
            
            try {
                // Algorithm parameters
                if (key == "populationSize") {
                    m_populationSize = std::stoi(value);
                    LOG_DEBUG("Setting population size to: " + std::to_string(m_populationSize));
                }
                else if (key == "mutationRate") {
                    m_mutationRate = std::stod(value);
                    LOG_DEBUG("Setting mutation rate to: " + std::to_string(m_mutationRate));
                }
                else if (key == "crossoverProbability") {
                    m_crossoverProbability = std::stod(value);
                    LOG_DEBUG("Setting crossover probability to: " + std::to_string(m_crossoverProbability));
                }
                else if (key == "targetFitness") {
                    m_targetFitness = std::stod(value);
                }
                else if (key == "tournamentSize") {
                    m_tournamentSize = std::stoi(value);
                    LOG_DEBUG("Setting tournament size to: " + std::to_string(m_tournamentSize));
                }
                else if (key == "replacementRatio") {
                    m_replacementRatio = std::stod(value);
                    LOG_DEBUG("Setting replacement ratio to: " + std::to_string(m_replacementRatio));
                }
                else if (key == "selectionMethod") {
                    m_selectionMethod = value;
                }
                else if (key == "noImprovementGenerations") {
                    m_noImprovementGenerations = std::stoi(value);
                    LOG_DEBUG("Setting no improvement generations to: " + std::to_string(m_noImprovementGenerations));
                }
                else if (key == "timeLimitSeconds") {
                    m_timeLimitSeconds = std::stoi(value);
                }
                // Instance-specific parameters
                else if (key == "k") {
                    m_k = std::stoi(value);
                    LOG_DEBUG("Setting k to: " + std::to_string(m_k));
                }
                else if (key == "deltaK") {
                    m_deltaK = std::stoi(value);
                    LOG_DEBUG("Setting deltaK to: " + std::to_string(m_deltaK));
                }
                else if (key == "lNeg") {
                    m_lNeg = std::stoi(value);
                    LOG_DEBUG("Setting lNeg to: " + std::to_string(m_lNeg));
                }
                else if (key == "lPoz") {
                    m_lPoz = std::stoi(value);
                    LOG_DEBUG("Setting lPoz to: " + std::to_string(m_lPoz));
                }
                else if (key == "repAllowed") {
                    m_repAllowed = (value == "true" || value == "1");
                    LOG_DEBUG("Setting repAllowed to: " + std::to_string(m_repAllowed));
                }
                else if (key == "probablePositive") {
                    m_probablePositive = std::stoi(value);
                    LOG_DEBUG("Setting probablePositive to: " + std::to_string(m_probablePositive));
                }
                // Adaptive parameters
                else if (key == "useAdaptiveMutation") {
                    m_adaptiveParams.useAdaptiveMutation = (value == "true");
                }
                else if (key == "minMutationRate") {
                    m_adaptiveParams.minMutationRate = std::stod(value);
                }
                else if (key == "maxMutationRate") {
                    m_adaptiveParams.maxMutationRate = std::stod(value);
                }
                else if (key == "stagnationGenerations") {
                    m_adaptiveParams.stagnationGenerations = std::stoi(value);
                }
                else if (key == "improvementThreshold") {
                    m_adaptiveParams.improvementThreshold = std::stod(value);
                }
                // Diversity parameters
                else if (key == "useFitnessSharing") {
                    m_diversityParams.useFitnessSharing = (value == "true");
                }
                else if (key == "useCrowding") {
                    m_diversityParams.useCrowding = (value == "true");
                }
                else if (key == "sharingRadius") {
                    m_diversityParams.sharingRadius = std::stod(value);
                }
                else if (key == "sharingAlpha") {
                    m_diversityParams.sharingAlpha = std::stod(value);
                }
                else if (key == "diversityWeight") {
                    m_diversityParams.diversityWeight = std::stod(value);
                }
                else if (key == "crossoverType") {
                    if (value == "adaptive") {
                        m_crossoverType = CrossoverType::ADAPTIVE;
                        LOG_DEBUG("Setting crossover type to adaptive");
                    }
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error parsing value for key " + key + ": " + e.what());
                return false;
            }
        }
    }
    
    LOG_INFO("Configuration loaded successfully from: " + filename);
    LOG_DEBUG("Final configuration values:");
    LOG_DEBUG("  Algorithm parameters:");
    LOG_DEBUG("    Population size: " + std::to_string(m_populationSize));
    LOG_DEBUG("    Mutation rate: " + std::to_string(m_mutationRate));
    LOG_DEBUG("    Crossover probability: " + std::to_string(m_crossoverProbability));
    LOG_DEBUG("    Tournament size: " + std::to_string(m_tournamentSize));
    LOG_DEBUG("    Replacement ratio: " + std::to_string(m_replacementRatio));
    LOG_DEBUG("    No improvement generations: " + std::to_string(m_noImprovementGenerations));
    LOG_DEBUG("  Instance parameters:");
    LOG_DEBUG("    k: " + std::to_string(m_k));
    LOG_DEBUG("    deltaK: " + std::to_string(m_deltaK));
    LOG_DEBUG("    lNeg: " + std::to_string(m_lNeg));
    LOG_DEBUG("    lPoz: " + std::to_string(m_lPoz));
    LOG_DEBUG("    repAllowed: " + std::to_string(m_repAllowed));
    LOG_DEBUG("    probablePositive: " + std::to_string(m_probablePositive));
    
    return true;
}

// ---------------------------------------
// Metody zwracające obiekty metaheurystyk
// ---------------------------------------
std::shared_ptr<IRepresentation> GAConfig::getRepresentation() const
{
    // Always return permutation representation
    return std::make_shared<PermutationRepresentation>();
}

std::shared_ptr<ISelection> GAConfig::getSelection() const
{
    return std::make_shared<TournamentSelection>(*this);
}

std::shared_ptr<ICrossover> GAConfig::getCrossover(const std::string& generationStr) const {
    if (m_crossoverType == CrossoverType::ADAPTIVE) {
        LOG_DEBUG("Getting adaptive crossover");
        if (!m_cachedAdaptiveCrossover) {
            LOG_DEBUG("Creating new adaptive crossover instance");
            m_cachedAdaptiveCrossover = std::make_shared<AdaptiveCrossover>(*this, m_instance);
        }
        // Set the generation if provided
        if (!generationStr.empty()) {
            try {
                int generation = std::stoi(generationStr);
                auto adaptiveCrossover = std::dynamic_pointer_cast<AdaptiveCrossover>(m_cachedAdaptiveCrossover);
                if (adaptiveCrossover) {
                    adaptiveCrossover->setGeneration(generation);
                    LOG_DEBUG("Set generation " + std::to_string(generation) + " for adaptive crossover");
                } else {
                    LOG_ERROR("Failed to cast to AdaptiveCrossover");
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Failed to parse generation number: " + generationStr);
            }
        }
        return m_cachedAdaptiveCrossover;
    }
    LOG_DEBUG("Using edge recombination crossover");
    return std::make_shared<EdgeRecombination>();
}

std::shared_ptr<IMutation> GAConfig::getMutation() const
{
    return std::make_shared<CombinedMutation>(m_mutationRate);
}

std::shared_ptr<IReplacement> GAConfig::getReplacement() const
{
    return std::make_shared<PartialReplacement>(m_replacementRatio, m_diversityParams.sharingRadius, m_cache);
}

std::shared_ptr<IFitness> GAConfig::getFitness() const {
    if (!m_cache) {
        throw std::runtime_error("Cache not set in GAConfig");
    }
    
    auto populationCache = std::dynamic_pointer_cast<PopulationCache>(m_cache);
    if (!populationCache) {
        throw std::runtime_error("Invalid cache type in GAConfig");
    }
    
    return std::make_shared<ConcreteOptimizedFitness>(*this, populationCache->getCurrentPopulation());
}

std::shared_ptr<IStopping> GAConfig::getStopping() const
{
    return std::make_shared<NoImprovementStopping>(m_noImprovementGenerations);
}

void GAConfig::setParameters(const ParameterSet& ps) {
    m_populationSize = ps.getInt("populationSize");
    m_mutationRate = ps.getDouble("mutationRate");
    m_crossoverProbability = ps.getDouble("crossoverRate");
    m_tournamentSize = ps.getInt("tournamentSize");
}

bool GAConfig::validate() const {
    std::shared_lock<std::shared_mutex> lock(configMutex);
    
    // Check numeric ranges
    if (m_populationSize <= 0) return false;
    if (m_mutationRate < 0.0 || m_mutationRate > 1.0) return false;
    if (m_replacementRatio < 0.0 || m_replacementRatio > 1.0) return false;
    if (m_crossoverProbability < 0.0 || m_crossoverProbability > 1.0) return false;
    if (m_tournamentSize <= 0) return false;
    if (m_noImprovementGenerations < 0) return false;
    if (m_timeLimitSeconds < 0) return false;
    
    // Check string parameters
    if (m_selectionMethod.empty()) return false;
    
    return true;
}

void GAConfig::setMutationRate(double rate) {
    // Clamp mutation rate to valid range for CombinedMutation
    static constexpr double MIN_MUTATION_RATE = 0.1;
    static constexpr double MAX_MUTATION_RATE = 0.4;
    
    if (rate < MIN_MUTATION_RATE || rate > MAX_MUTATION_RATE) {
        LOG_WARNING("Mutation rate " + std::to_string(rate) + 
                   " outside valid range [" + std::to_string(MIN_MUTATION_RATE) + 
                   ", " + std::to_string(MAX_MUTATION_RATE) + "] - clamping");
    }
    m_mutationRate = std::clamp(rate, MIN_MUTATION_RATE, MAX_MUTATION_RATE);
}

