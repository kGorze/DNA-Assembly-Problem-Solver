#pragma once

#include <memory>
#include <string>

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

    // -------------------------
    // Getters for required fields
    // -------------------------
    int getPopulationSize() const { return populationSize; }
    int getMaxGenerations()  const { return maxGenerations; }

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
