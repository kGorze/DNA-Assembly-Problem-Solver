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
 * Adaptive crossover parameters structure
 */
struct AdaptiveCrossoverParams {
    double inertia{0.7};
    int adaptationInterval{20};
    int minTrials{5};
    double minProb{0.1};
};

class GAConfig {
public:
    GAConfig() = default;
    GAConfig(const GAConfig& other) = default;
    GAConfig& operator=(const GAConfig& other) = default;

    // Load configuration from file
    bool loadFromFile(const std::string& filename);
    
    // Reset to default values
    void resetToDefaults();

    // Getters
    int getPopulationSize() const { return m_populationSize; }
    double getMutationRate() const { return m_mutationRate; }
    double getCrossoverProbability() const { return m_crossoverProbability; }
    double getTargetFitness() const { return m_targetFitness; }
    int getTournamentSize() const { return m_tournamentSize; }
    int getK() const { return m_k; }
    int getDeltaK() const { return m_deltaK; }
    int getLNeg() const { return m_lNeg; }
    int getLPoz() const { return m_lPoz; }
    bool isRepAllowed() const { return m_repAllowed; }
    double getProbablePositive() const { return m_probablePositive; }
    double getReplacementRatio() const { return m_replacementRatio; }
    std::string getSelectionMethod() const { return m_selectionMethod; }
    int getNoImprovementGenerations() const { return m_noImprovementGenerations; }
    int getTimeLimitSeconds() const { return m_timeLimitSeconds; }

    // Setters
    void setPopulationSize(int size) { m_populationSize = size; }
    void setMutationRate(double rate) { m_mutationRate = rate; }
    void setCrossoverProbability(double probability) { m_crossoverProbability = probability; }
    void setTargetFitness(double fitness) { m_targetFitness = fitness; }
    void setTournamentSize(int size) { m_tournamentSize = size; }
    void setK(int k) { m_k = k; }
    void setDeltaK(int deltaK) { m_deltaK = deltaK; }
    void setLNeg(int lNeg) { m_lNeg = lNeg; }
    void setLPoz(int lPoz) { m_lPoz = lPoz; }
    void setRepAllowed(bool allowed) { m_repAllowed = allowed; }
    void setProbablePositive(double prob) { m_probablePositive = prob; }
    void setReplacementRatio(double ratio) { m_replacementRatio = ratio; }
    void setSelectionMethod(const std::string& method) { m_selectionMethod = method; }
    void setNoImprovementGenerations(int gens) { m_noImprovementGenerations = gens; }
    void setTimeLimitSeconds(int seconds) { m_timeLimitSeconds = seconds; }

    // Cache management
    void setCache(std::shared_ptr<IPopulationCache> cachePtr) {
        m_cache = cachePtr;
    }
    
    std::shared_ptr<IPopulationCache> getCache() const {
        return m_cache;
    }

    // Component getters - these always return the same implementations
    std::shared_ptr<IRepresentation> getRepresentation() const;  // Always returns PermutationRepresentation
    std::shared_ptr<ISelection> getSelection() const;
    std::shared_ptr<ICrossover> getCrossover(const std::string& = "") const;  // Always returns AdaptiveCrossover
    std::shared_ptr<IMutation> getMutation() const;  // Always returns OnePointMutation
    std::shared_ptr<IReplacement> getReplacement() const;
    std::shared_ptr<IFitness> getFitness() const;  // Always returns OptimizedGraphBasedFitness
    std::shared_ptr<IStopping> getStopping() const;  // Always returns NoImprovementStopping with time limit

    // Parameter management
    void setParameters(const ParameterSet& ps);
    bool validate() const;

    const AdaptiveCrossoverParams& getAdaptiveParams() const {
        return adaptiveParams;
    }

    int getParentCount() const { return 2; }

private:
    int m_populationSize = 100;
    double m_mutationRate = 0.1;
    double m_crossoverProbability = 0.8;
    double m_targetFitness = 1.0;
    int m_tournamentSize = 5;
    int m_k = 0;
    int m_deltaK = 0;
    int m_lNeg = 0;
    int m_lPoz = 0;
    bool m_repAllowed = false;
    double m_probablePositive = 0.0;
    double m_replacementRatio = 0.5;
    std::string m_selectionMethod = "tournament";
    int m_noImprovementGenerations = 30;
    int m_timeLimitSeconds = 60;
    
    AdaptiveCrossoverParams adaptiveParams;
    
    // Cache for population
    std::shared_ptr<IPopulationCache> m_cache;
};
