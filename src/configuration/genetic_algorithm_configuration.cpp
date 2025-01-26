//
// Created by konrad_guest on 23/01/2025.
//
#include "../../include/configuration/genetic_algorithm_configuration.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
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
    m_mutationRate = 0.1;
    m_crossoverProbability = 0.8;
    m_targetFitness = 1.0;
    m_tournamentSize = 5;
    m_replacementRatio = 0.5;
    m_selectionMethod = "tournament";
    m_k = 8;
    m_deltaK = 2;
    m_lNeg = 25;
    m_lPoz = 25;
    m_repAllowed = false;
    m_probablePositive = 0;
    
    // Initialize adaptive crossover parameters
    adaptiveParams.inertia = 0.7;
    adaptiveParams.adaptationInterval = 20;
    adaptiveParams.minTrials = 5;
    adaptiveParams.minProb = 0.1;
    
    // Set stopping criteria parameters
    m_noImprovementGenerations = 30;  // Default to 30 generations without improvement
    m_timeLimitSeconds = 60;          // Default to 60 seconds time limit
}

bool GAConfig::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        LOG_ERROR("Failed to open configuration file: " + filename);
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        std::string value;

        if (std::getline(iss, key, '=') && std::getline(iss, value)) {
            // Trim whitespace from key and value
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);
            
            try {
                if (key == "populationSize") m_populationSize = std::stoi(value);
                else if (key == "mutationRate") m_mutationRate = std::stod(value);
                else if (key == "crossoverProbability") m_crossoverProbability = std::stod(value);
                else if (key == "targetFitness") m_targetFitness = std::stod(value);
                else if (key == "tournamentSize") m_tournamentSize = std::stoi(value);
                else if (key == "replacementRatio") m_replacementRatio = std::stod(value);
                else if (key == "selectionMethod") m_selectionMethod = value;
                else if (key == "k") m_k = std::stoi(value);
                else if (key == "deltaK") m_deltaK = std::stoi(value);
                else if (key == "lNeg") m_lNeg = std::stoi(value);
                else if (key == "lPoz") m_lPoz = std::stoi(value);
                else if (key == "repAllowed") m_repAllowed = (value == "true" || value == "1");
                else if (key == "probablePositive") m_probablePositive = std::stod(value);
                else if (key == "noImprovementGenerations") m_noImprovementGenerations = std::stoi(value);
                else if (key == "timeLimitSeconds") m_timeLimitSeconds = std::stoi(value);
                // Adaptive crossover parameters
                else if (key == "adaptiveInertia") adaptiveParams.inertia = std::stod(value);
                else if (key == "adaptationInterval") adaptiveParams.adaptationInterval = std::stoi(value);
                else if (key == "minTrials") adaptiveParams.minTrials = std::stoi(value);
                else if (key == "minProb") adaptiveParams.minProb = std::stod(value);
                else {
                    LOG_WARNING("Unknown configuration key: " + key);
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error parsing value for key " + key + ": " + e.what());
                return false;
            }
        }
    }

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

std::shared_ptr<ICrossover> GAConfig::getCrossover(const std::string&) const
{
    // Always return adaptive crossover
    return std::make_shared<AdaptiveCrossover>(*this);
}

std::shared_ptr<IMutation> GAConfig::getMutation() const
{
    return std::make_shared<PointMutation>(m_mutationRate);
}

std::shared_ptr<IReplacement> GAConfig::getReplacement() const
{
    return std::make_shared<PartialReplacement>(m_replacementRatio, m_cache);
}

std::shared_ptr<IFitness> GAConfig::getFitness() const
{
    // Always return optimized graph-based fitness
    return std::make_shared<OptimizedGraphBasedFitness>();
}

std::shared_ptr<IStopping> GAConfig::getStopping() const
{
    return std::make_shared<NoImprovementStopping>(50);  // Default to 50 generations without improvement
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

