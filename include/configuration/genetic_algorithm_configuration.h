#pragma once

#include <memory>
#include <string>
#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"
#include "metaheuristics/population_cache.h"
#include "metaheuristics/representation.h"
#include "metaheuristics/adaptive_crossover.h"



// Singleton GAConfig class
class GAConfig {
public:
    // Singleton access
    static GAConfig& getInstance();

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
    GAConfig(); // Private constructor for singleton

    std::shared_ptr<IPopulationCache> m_cache;
};
