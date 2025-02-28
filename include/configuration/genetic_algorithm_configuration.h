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
 * Adaptive parameters structure
 */
struct AdaptiveParams {
    bool useAdaptiveMutation = true;
    double minMutationRate = 0.1;
    double maxMutationRate = 0.4;
    int stagnationGenerations = 5;  // Generations without improvement before adaptation
    double improvementThreshold = 0.01; // Minimum improvement to reset stagnation counter
};

/**
 * Diversity parameters structure
 */
struct DiversityParams {
    bool useFitnessSharing = true;
    bool useCrowding = false;
    double sharingRadius = 0.2;  // Radius for fitness sharing
    double sharingAlpha = 1.0;   // Shape parameter for sharing function
    double diversityWeight = 0.3; // Weight for diversity score in fitness
};

// Add before GAConfig class
enum class CrossoverType {
    ADAPTIVE,
    // Add other types if needed
};

class GAConfig {
public:
    GAConfig();
    virtual ~GAConfig() = default;

    // Load configuration from file
    bool loadFromFile(const std::string& filename);
    
    // Reset to default values
    void resetToDefaults();

    // Getters
    int getPopulationSize() const { return m_populationSize; }
    int getMaxGenerations() const { return m_maxGenerations; }
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
    const std::string& getSelectionMethod() const { return m_selectionMethod; }
    int getNoImprovementGenerations() const { return m_noImprovementGenerations; }
    int getTimeLimitSeconds() const { return m_timeLimitSeconds; }
    const AdaptiveParams& getAdaptiveParams() const { return m_adaptiveParams; }
    const DiversityParams& getDiversityParams() const { return m_diversityParams; }
    CrossoverType getCrossoverType() const { return m_crossoverType; }
    const std::shared_ptr<DNAInstance>& getInstance() const { return m_instance; }

    // Setters
    void setPopulationSize(int size) { m_populationSize = size; }
    void setMaxGenerations(int generations) { m_maxGenerations = generations; }
    void setMutationRate(double rate);  // Defined in cpp file
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
    void setAdaptiveParams(const AdaptiveParams& params) { m_adaptiveParams = params; }
    void setDiversityParams(const DiversityParams& params) { m_diversityParams = params; }
    void setCrossoverType(CrossoverType type) { m_crossoverType = type; }
    void setInstance(const std::shared_ptr<DNAInstance>& instance) { m_instance = instance; }

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

    int getParentCount() const { return 2; }

protected:
    AdaptiveParams m_adaptiveParams;
    DiversityParams m_diversityParams;

private:
    int m_populationSize = 100;
    int m_maxGenerations = 1000;
    double m_mutationRate = 0.1;
    double m_crossoverProbability = 0.8;
    double m_targetFitness = 1.0;
    int m_tournamentSize = 2;
    int m_k = 7;
    int m_deltaK = 0;
    int m_lNeg = 10;
    int m_lPoz = 10;
    bool m_repAllowed = true;
    double m_probablePositive = 0.0;
    double m_replacementRatio = 0.7;
    std::string m_selectionMethod = "rank";
    int m_noImprovementGenerations = 30;
    int m_timeLimitSeconds = 60;
    
    // Cache for population
    std::shared_ptr<IPopulationCache> m_cache;
    CrossoverType m_crossoverType;
    std::shared_ptr<DNAInstance> m_instance;
    mutable std::shared_ptr<ICrossover> m_cachedAdaptiveCrossover;
};
