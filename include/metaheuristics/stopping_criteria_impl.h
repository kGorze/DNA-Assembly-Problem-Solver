#pragma once

#include "../interfaces/i_stopping.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "../utils/logging.h"
#include <vector>
#include <memory>
#include <limits>
#include <mutex>
#include <chrono>

class NoImprovementStopping : public IStopping {
public:
    NoImprovementStopping(int maxGenerationsWithoutImprovement = 50)
        : m_maxGenerationsWithoutImprovement(maxGenerationsWithoutImprovement)
        , m_generationsWithoutImprovement(0)
        , m_bestFitnessSoFar(-1e9) {}

    bool shouldStop([[maybe_unused]] int currentGeneration, double bestFitness) const override {
        if (bestFitness > m_bestFitnessSoFar) {
            m_bestFitnessSoFar = bestFitness;
            m_generationsWithoutImprovement = 0;
        } else {
            m_generationsWithoutImprovement++;
        }
        return m_generationsWithoutImprovement >= m_maxGenerationsWithoutImprovement;
    }

    void reset() override {
        m_bestFitnessSoFar = -1e9;
        m_generationsWithoutImprovement = 0;
        LOG_INFO("NoImprovementStopping criteria reset");
    }

private:
    const int m_maxGenerationsWithoutImprovement;
    mutable int m_generationsWithoutImprovement;
    mutable double m_bestFitnessSoFar;
};

class MaxGenerationsStopping : public IStopping {
public:
    MaxGenerationsStopping(int maxGenerations = 100)
        : m_maxGenerations(maxGenerations) {}

    bool shouldStop(int currentGeneration, [[maybe_unused]] double bestFitness) const override {
        return currentGeneration >= m_maxGenerations;
    }

    void reset() override {
        // Nothing to reset for MaxGenerationsStopping
    }

private:
    const int m_maxGenerations;
};

class TimeLimitStopping : public IStopping {
public:
    explicit TimeLimitStopping(int timeLimitInSeconds = 60) 
        : m_timeLimit(timeLimitInSeconds) {
        m_startTime = std::chrono::high_resolution_clock::now();
    }

    bool shouldStop([[maybe_unused]] int currentGeneration, 
                   [[maybe_unused]] double bestFitness) const override {
        auto currentTime = std::chrono::high_resolution_clock::now();
        auto elapsedTime = std::chrono::duration_cast<std::chrono::seconds>(currentTime - m_startTime);
        return elapsedTime >= m_timeLimit;
    }

    void reset() override {
        m_startTime = std::chrono::high_resolution_clock::now();
    }

private:
    std::chrono::seconds m_timeLimit;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_startTime;
}; 