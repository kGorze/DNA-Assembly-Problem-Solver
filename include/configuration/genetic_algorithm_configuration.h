#pragma once

#include <memory>
#include <string>
#include <iostream>

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
    // Singleton access
    static GAConfig& getInstance();

    double getGlobalBestFitness() const;
    void setGlobalBestFitness(double fitness);

    // -------------------------
    // Public parameters
    // -------------------------
    int populationSize;      
    double mutationRate;     
    double replacementRatio;
    int maxGenerations;      

    double crossoverProbability;

    std::string selectionMethod;
    std::string crossoverType;
    std::string mutationMethod;
    std::string replacementMethod;
    std::string stoppingMethod;
    int noImprovementGenerations;
    int tournamentSize;
    int timeLimitSeconds;

    struct AdaptiveCrossoverParams {
        double inertia;
        int adaptationInterval;
        int minTrials;
        double minProb;
    } adaptiveParams;

    std::string fitnessType;
    double alpha;
    double beta;

    // DNA Generation parameters
    int k;                  // długość oligo
    int deltaK;            // maksymalna zmiana długości oligo
    int lNeg;              // liczba błędów negatywnych
    int lPoz;              // liczba błędów pozytywnych
    bool repAllowed;       // czy dozwolone powtórzenia w DNA
    int probablePositive;  // sposób generowania błędów pozytywnych

    // -------------------------
    // Getters for required fields
    // -------------------------
    int getPopulationSize() const { return populationSize; }
    int getMaxGenerations()  const { return maxGenerations; }
    int getTournamentSize() const { return tournamentSize; }

    void setParameters(const ParameterSet &ps) {
        // sprawdzaj klucze i konwertuj
        // (ważne, by użyć try/catch lub sprawdzać, czy klucz istnieje)
        if (ps.params.count("populationSize")) {
            populationSize = std::stoi(ps.params.at("populationSize"));
        }
        if (ps.params.count("mutationRate")) {
            mutationRate = std::stod(ps.params.at("mutationRate"));
        }
        if (ps.params.count("replacementRatio")) {
            replacementRatio = std::stod(ps.params.at("replacementRatio"));
        }
        if (ps.params.count("tournamentSize")) {
            tournamentSize = std::stoi(ps.params.at("tournamentSize"));
        }
        if (ps.params.count("crossoverType")) {
            crossoverType = ps.params.at("crossoverType");
        }
        if (ps.params.count("selectionMethod")) {
            selectionMethod = ps.params.at("selectionMethod");
        }
        if (ps.params.count("adaptive.inertia")) {
            adaptiveParams.inertia = std::stod(ps.params.at("adaptive.inertia"));
        }
        if (ps.params.count("adaptive.adaptationInterval")) {
            adaptiveParams.adaptationInterval = std::stoi(ps.params.at("adaptive.adaptationInterval"));
        }
        if (ps.params.count("adaptive.minTrials")) {
            adaptiveParams.minTrials = std::stoi(ps.params.at("adaptive.minTrials"));
        }
        if (ps.params.count("adaptive.minProb")) {
            adaptiveParams.minProb = std::stod(ps.params.at("adaptive.minProb"));
        }
        if (ps.params.count("k")) {
            k = std::stoi(ps.params.at("k"));
        }
        if (ps.params.count("deltaK")) {
            deltaK = std::stoi(ps.params.at("deltaK"));
        }
        if (ps.params.count("lNeg")) {
            lNeg = std::stoi(ps.params.at("lNeg"));
        }
        if (ps.params.count("lPoz")) {
            lPoz = std::stoi(ps.params.at("lPoz"));
        }
        if (ps.params.count("repAllowed")) {
            repAllowed = ps.params.at("repAllowed") == "true";
        }
        if (ps.params.count("probablePositive")) {
            probablePositive = std::stoi(ps.params.at("probablePositive"));
        }
    }


    // -------------------------
    // Interface to create operator objects
    // -------------------------
    std::shared_ptr<IRepresentation> getRepresentation() const;
    std::shared_ptr<ISelection>      getSelection()      const;
    std::shared_ptr<ICrossover>      getCrossover(const std::string& type) const;
    std::shared_ptr<IMutation>       getMutation()       const;
    std::shared_ptr<IReplacement>    getReplacement()    const;
    std::shared_ptr<IFitness>        getFitness()        const;
    std::shared_ptr<IStopping>       getStopping()       const;

    // -------------------------
    // Cache management
    // -------------------------
    void setCache(std::shared_ptr<IPopulationCache> cache) { m_cache = cache; }
    std::shared_ptr<IPopulationCache> getCache() const { return m_cache; }

    // -------------------------
    // Configuration loading
    // -------------------------
    bool loadFromFile(const std::string& filePath);

    // -------------------------
    // Setters with validation
    // -------------------------
    void setReplacementRatio(double ratio);

    bool validate() const {
        if (k <= 0) {
            std::cerr << "[ERROR] K must be positive" << std::endl;
            return false;
        }
        
        if (populationSize <= 0) {
            std::cerr << "[ERROR] Population size must be positive" << std::endl;
            return false;
        }
        
        if (mutationRate < 0.0 || mutationRate > 1.0) {
            std::cerr << "[ERROR] Mutation rate must be between 0 and 1" << std::endl;
            return false;
        }
        
        if (replacementRatio < 0.0 || replacementRatio > 1.0) {
            std::cerr << "[ERROR] Replacement ratio must be between 0 and 1" << std::endl;
            return false;
        }
        
        return true;
    }

private:
    GAConfig();
    GAConfig(const GAConfig&) = delete;
    GAConfig& operator=(const GAConfig&) = delete;

    // Zmienna, w której będziemy przechowywać najlepszy fitness
    double m_globalBestFit = -std::numeric_limits<double>::infinity();
    double globalBestFitness;
    // Cache do populacji
    std::shared_ptr<IPopulationCache> m_cache;
};
