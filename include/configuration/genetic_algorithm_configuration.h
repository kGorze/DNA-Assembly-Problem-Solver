#pragma once

#include <memory>
#include <string>
#include <iostream>
#include <mutex>

#include <limits>
#include <unordered_map>
#include "metaheuristics/representation.h"
#include "../include/tuning/tuning_structures.h"

#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"
#include "metaheuristics/population_cache.h"
#include "metaheuristics/representation.h"
#include "metaheuristics/adaptive_crossover.h"

// Add these forward declarations at the top of the file
class ISelection;
class ICrossover;
class IMutation;
class IReplacement;
class IFitness;
class IStopping;
class IRepresentation;
class IPopulationCache;

/**

 * Struktura parametrów adaptacyjnego krzyżowania.

 */

struct AdaptiveCrossoverParams {

    double inertia;
    int adaptationInterval;
    int minTrials;
    double minProb;
};


// Singleton GAConfig class
class GAConfig {
public:
    // Make constructor public
    GAConfig();
    
    // Configuration loading
    bool loadFromFile(const std::string& filePath);
    void resetToDefaults();
    
    // Getters for components
    std::shared_ptr<IRepresentation> getRepresentation() const;
    std::shared_ptr<ISelection> getSelection() const;
    std::shared_ptr<ICrossover> getCrossover(const std::string& type) const;
    std::shared_ptr<IMutation> getMutation() const;
    std::shared_ptr<IReplacement> getReplacement() const;
    std::shared_ptr<IFitness> getFitness() const;
    std::shared_ptr<IStopping> getStopping() const;
    
    // Getters/setters for parameters
    void setMaxGenerations(int value) { 
        std::lock_guard<std::mutex> lock(configMutex);
        m_maxGenerations = value;
        std::cout << "[GAConfig] Setting maxGenerations = " << value << std::endl;
    }
    int getMaxGenerations() const { 
        std::lock_guard<std::mutex> lock(configMutex);
        std::cout << "[GAConfig] Getting maxGenerations = " << m_maxGenerations << std::endl;
        return m_maxGenerations;
    }
    void setReplacementRatio(double ratio);
    
    // Global best fitness tracking
    double getGlobalBestFitness() const;
    void setGlobalBestFitness(double fitness);
    
    // Parameter validation
    bool validate() const;
    
    // Cache management
    void setCache(std::shared_ptr<IPopulationCache> cache) { m_cache = cache; }
    
    // Parameter setting from ParameterSet
    void setParameters(const ParameterSet &ps);

    // Public member variables (consider making these private with getters/setters)
    int populationSize;
    double mutationRate;
    double replacementRatio;
    double crossoverProbability;
    
    std::string selectionMethod;
    std::string crossoverType;
    std::string mutationMethod;
    std::string replacementMethod;
    std::string stoppingMethod;
    std::string fitnessType;
    
    int noImprovementGenerations;
    int tournamentSize;
    int timeLimitSeconds;
    
    AdaptiveCrossoverParams adaptiveParams;
    
    double alpha;
    double beta;
    
    // DNA Generation parameters
    int k;
    int deltaK;
    int lNeg;
    int lPoz;
    bool repAllowed;
    int probablePositive;

private:
    mutable std::mutex configMutex;
    int m_maxGenerations;
    double m_globalBestFit;
    std::string lastLoadedConfig;
    bool isInitialized;
    std::shared_ptr<IPopulationCache> m_cache;
};
