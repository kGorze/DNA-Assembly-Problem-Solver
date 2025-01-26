#pragma once

#include <memory>
#include <string>
#include <iostream>
#include <mutex>
#include <atomic>
#include <shared_mutex>
#include <limits>
#include <unordered_map>

// Interface includes
#include "../interfaces/i_representation.h"
#include "../interfaces/i_selection.h"
#include "../interfaces/i_crossover.h"
#include "../interfaces/i_mutation.h"
#include "../interfaces/i_replacement.h"
#include "../interfaces/i_fitness.h"
#include "../interfaces/i_stopping.h"
#include "../interfaces/i_population_cache.h"

// Other includes
#include "../tuning/tuning_structures.h"

// Forward declarations
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
    double inertia{0.7};
    int adaptationInterval{20};
    int minTrials{5};
    double minProb{0.1};
};

class GAConfig {
public:
    // Make constructor public and remove singleton
    GAConfig() { resetToDefaults(); }
    
    // Load configuration from file
    bool loadFromFile(const std::string& filename);
    
    // Reset to default values
    void resetToDefaults();

    // Getters with thread safety
    int getPopulationSize() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return m_populationSize;
    }

    double getMutationRate() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return m_mutationRate;
    }

    // Setters with thread safety and validation
    void setPopulationSize(int size) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (size > 0) {
            m_populationSize = size;
        }
    }

    void setMutationRate(double rate) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (rate >= 0.0 && rate <= 1.0) {
            m_mutationRate = rate;
        }
    }
    
    // Thread-safe getters/setters using shared mutex for better read performance
    void setMaxGenerations(int value) { 
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (value <= 0) {
            std::cerr << "[GAConfig] Invalid maxGenerations value: " << value << ", using default of 100" << std::endl;
            m_maxGenerations = 100;
        } else {
            m_maxGenerations = value;
        }
        std::cout << "[GAConfig] Setting maxGenerations = " << m_maxGenerations << std::endl;
    }
    
    int getMaxGenerations() const { 
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return m_maxGenerations;
    }
    
    double getReplacementRatio() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return replacementRatio;
    }
    
    double getCrossoverProbability() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return crossoverProbability;
    }
    
    std::string getSelectionMethod() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return selectionMethod;
    }
    
    std::string getCrossoverType() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return crossoverType;
    }
    
    std::string getMutationMethod() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return mutationMethod;
    }

    // Cache management
    void setCache(std::shared_ptr<IPopulationCache> cachePtr) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        m_cache = cachePtr;
    }
    
    std::shared_ptr<IPopulationCache> getCache() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return m_cache;
    }

    // Component getters
    std::shared_ptr<IRepresentation> getRepresentation() const;
    std::shared_ptr<ISelection> getSelection() const;
    std::shared_ptr<ICrossover> getCrossover(const std::string& type) const;
    std::shared_ptr<IMutation> getMutation() const;
    std::shared_ptr<IReplacement> getReplacement() const;
    std::shared_ptr<IFitness> getFitness() const;
    std::shared_ptr<IStopping> getStopping() const;

    // Instance parameters
    void setK(int value) { k = value; }
    void setDeltaK(int value) { deltaK = value; }
    void setLNeg(int value) { lNeg = value; }
    void setLPoz(int value) { lPoz = value; }
    void setRepAllowed(bool value) { repAllowed = value; }
    void setProbablePositive(int value) { probablePositive = value; }
    
    int getK() const { return k; }
    int getDeltaK() const { return deltaK; }
    int getLNeg() const { return lNeg; }
    int getLPoz() const { return lPoz; }
    bool isRepAllowed() const { return repAllowed; }
    int getProbablePositive() const { return probablePositive; }

    // Add missing getters and setters
    void setCrossoverProbability(double value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (value >= 0.0 && value <= 1.0) crossoverProbability = value;
    }

    void setReplacementRatio(double ratio) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (ratio >= 0.0 && ratio <= 1.0) replacementRatio = ratio;
    }

    int getTournamentSize() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return tournamentSize;
    }

    void setTournamentSize(int size) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (size > 0) tournamentSize = size;
    }

    // Global best fitness tracking
    double getGlobalBestFitness() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return m_globalBestFit;
    }

    void setGlobalBestFitness(double fitness) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        m_globalBestFit = fitness;
    }

    // Parameter management
    void setParameters(const ParameterSet& ps);
    bool validate() const;

    const AdaptiveCrossoverParams& getAdaptiveParams() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return adaptiveParams;
    }

    int getParentCount() const { return 2; }
    double getTargetFitness() const { return targetFitness; }

private:
    mutable std::shared_mutex configMutex;
    std::atomic<bool> isInitialized{false};
    std::string lastLoadedConfig;
    
    // Configuration parameters
    int m_maxGenerations{100};
    int m_populationSize{100};
    double m_mutationRate{0.15};
    double replacementRatio{0.7};
    double crossoverProbability{1.0};
    
    std::string selectionMethod{"tournament"};
    std::string crossoverType{"order"};
    std::string mutationMethod{"point"};
    std::string replacementMethod{"partial"};
    std::string stoppingMethod{"maxGenerations"};
    std::string fitnessType{"optimized_graph"};
    
    int noImprovementGenerations{30};
    int tournamentSize{3};
    int timeLimitSeconds{60};
    
    AdaptiveCrossoverParams adaptiveParams;
    
    double alpha{0.7};
    double beta{0.3};
    
    // DNA Generation parameters
    int k{8};
    int deltaK{2};
    int lNeg{25};
    int lPoz{25};
    bool repAllowed{false};
    int probablePositive{0};
    
    // Global best fitness tracking
    double m_globalBestFit{-std::numeric_limits<double>::infinity()};
    
    // Cache for population
    std::shared_ptr<IPopulationCache> m_cache;

    double targetFitness;
};
