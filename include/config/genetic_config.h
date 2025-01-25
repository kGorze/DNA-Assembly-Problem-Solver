#pragma once

#include <string>

class GeneticConfig {
public:
    // Constructor with default values
    GeneticConfig(
        int populationSize = 100,
        double mutationRate = 0.1,
        double crossoverRate = 0.8,
        double replacementRatio = 0.7,
        int tournamentSize = 3,
        int maxGenerations = 1000,
        int logInterval = 10,
        int cacheCleanupInterval = 100
    ) : m_populationSize(populationSize)
      , m_mutationRate(mutationRate)
      , m_crossoverRate(crossoverRate)
      , m_replacementRatio(replacementRatio)
      , m_tournamentSize(tournamentSize)
      , m_maxGenerations(maxGenerations)
      , m_logInterval(logInterval)
      , m_cacheCleanupInterval(cacheCleanupInterval)
    {}

    // Getters
    int getPopulationSize() const { return m_populationSize; }
    double getMutationRate() const { return m_mutationRate; }
    double getCrossoverRate() const { return m_crossoverRate; }
    double getReplacementRatio() const { return m_replacementRatio; }
    int getTournamentSize() const { return m_tournamentSize; }
    int getMaxGenerations() const { return m_maxGenerations; }
    int getLogInterval() const { return m_logInterval; }
    int getCacheCleanupInterval() const { return m_cacheCleanupInterval; }
    double getGlobalBestFitness() const { return m_globalBestFitness; }

    // Setters
    void setPopulationSize(int value) { m_populationSize = value; }
    void setMutationRate(double value) { m_mutationRate = value; }
    void setCrossoverRate(double value) { m_crossoverRate = value; }
    void setReplacementRatio(double value) { m_replacementRatio = value; }
    void setTournamentSize(int value) { m_tournamentSize = value; }
    void setMaxGenerations(int value) { m_maxGenerations = value; }
    void setLogInterval(int value) { m_logInterval = value; }
    void setCacheCleanupInterval(int value) { m_cacheCleanupInterval = value; }
    void setGlobalBestFitness(double value) { m_globalBestFitness = value; }

private:
    int m_populationSize;
    double m_mutationRate;
    double m_crossoverRate;
    double m_replacementRatio;
    int m_tournamentSize;
    int m_maxGenerations;
    int m_logInterval;
    int m_cacheCleanupInterval;
    double m_globalBestFitness = -std::numeric_limits<double>::infinity();
}; 