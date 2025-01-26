#pragma once

class GAConfig {
public:
    GAConfig() = default;

    // Getters
    int getPopulationSize() const { return m_populationSize; }
    int getMaxGenerations() const { return m_maxGenerations; }
    double getMutationRate() const { return m_mutationRate; }
    double getCrossoverProbability() const { return m_crossoverProbability; }
    double getTargetFitness() const { return m_targetFitness; }
    int getTournamentSize() const { return m_tournamentSize; }

    // Setters
    void setPopulationSize(int size) { m_populationSize = size; }
    void setMaxGenerations(int generations) { m_maxGenerations = generations; }
    void setMutationRate(double rate) { m_mutationRate = rate; }
    void setCrossoverProbability(double probability) { m_crossoverProbability = probability; }
    void setTargetFitness(double fitness) { m_targetFitness = fitness; }
    void setTournamentSize(int size) { m_tournamentSize = size; }

private:
    int m_populationSize = 100;
    int m_maxGenerations = 1000;
    double m_mutationRate = 0.1;
    double m_crossoverProbability = 0.8;
    double m_targetFitness = 1.0;
    int m_tournamentSize = 5;
}; 