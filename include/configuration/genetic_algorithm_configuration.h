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

    double inertia;
    int adaptationInterval;
    int minTrials;
    double minProb;
};


// Singleton GAConfig class
class GAConfig {
public:
    // Make constructor public and return singleton instance
    static GAConfig& getInstance() {
        static GAConfig instance;
        return instance;
    }

    // Load configuration from file
    bool loadFromFile(const std::string& filename);

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

    // Configuration loading - now returns void since we'll use exceptions for errors
    void resetToDefaults();
    
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
    
    std::string getReplacementMethod() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return replacementMethod;
    }
    
    std::string getStoppingMethod() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return stoppingMethod;
    }
    
    std::string getFitnessType() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return fitnessType;
    }
    
    int getNoImprovementGenerations() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return noImprovementGenerations;
    }
    
    int getTournamentSize() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return tournamentSize;
    }
    
    int getTimeLimitSeconds() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return timeLimitSeconds;
    }
    
    AdaptiveCrossoverParams getAdaptiveParams() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return adaptiveParams;
    }
    
    double getAlpha() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return alpha;
    }
    
    double getBeta() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return beta;
    }
    
    // DNA Generation parameters getters
    int getK() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return k;
    }
    
    int getDeltaK() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return deltaK;
    }
    
    int getLNeg() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return lNeg;
    }
    
    int getLPoz() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return lPoz;
    }
    
    bool getRepAllowed() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return repAllowed;
    }
    
    int getProbablePositive() const {
        std::shared_lock<std::shared_mutex> lock(configMutex);
        return probablePositive;
    }
    
    // Setters with validation
    void setK(int value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (value > 0) k = value;
    }
    
    void setDeltaK(int value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        deltaK = value;
    }
    
    void setLNeg(int value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (value >= 0) lNeg = value;
    }
    
    void setLPoz(int value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (value >= 0) lPoz = value;
    }
    
    void setRepAllowed(bool value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        repAllowed = value;
    }
    
    void setProbablePositive(int value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        probablePositive = value;
    }
    
    void setCrossoverProbability(double value) {
        std::unique_lock<std::shared_mutex> lock(configMutex);
        if (value >= 0.0 && value <= 1.0) crossoverProbability = value;
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

    // Keep existing component getters
    std::shared_ptr<IRepresentation> getRepresentation() const;
    std::shared_ptr<ISelection> getSelection() const;
    std::shared_ptr<ICrossover> getCrossover(const std::string& type) const;
    std::shared_ptr<IMutation> getMutation() const;
    std::shared_ptr<IReplacement> getReplacement() const;
    std::shared_ptr<IFitness> getFitness() const;
    std::shared_ptr<IStopping> getStopping() const;

    // Add missing method declarations
    void setReplacementRatio(double ratio);
    double getGlobalBestFitness() const;
    void setGlobalBestFitness(double fitness);
    void setParameters(const ParameterSet& ps);
    bool validate() const;

private:
    GAConfig() = default;
    GAConfig(const GAConfig&) = delete;
    GAConfig& operator=(const GAConfig&) = delete;
    
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
    
    AdaptiveCrossoverParams adaptiveParams{0.7, 20, 5, 0.1};
    
    double alpha{0.7};
    double beta{0.3};
    
    // DNA Generation parameters
    int k{8};
    int deltaK{2};
    int lNeg{25};
    int lPoz{25};
    bool repAllowed{false};
    int probablePositive{0};
    
    double m_globalBestFit{-std::numeric_limits<double>::infinity()};
    std::shared_ptr<IPopulationCache> m_cache;
};
