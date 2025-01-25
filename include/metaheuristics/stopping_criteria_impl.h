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
    explicit NoImprovementStopping(int maxGenerationsWithoutImprovement)
        : m_maxGenerationsWithoutImprovement(std::max(1, maxGenerationsWithoutImprovement))
        , m_bestFitness(-std::numeric_limits<double>::infinity())
        , m_generationsWithoutImprovement(0)
    {
        if (maxGenerationsWithoutImprovement <= 0) {
            LOG_WARNING("Invalid maxGenerationsWithoutImprovement value: " + 
                       std::to_string(maxGenerationsWithoutImprovement) + 
                       ". Using default value of 50");
            m_maxGenerationsWithoutImprovement = 50;
        }
    }

    bool shouldStop(int currentGeneration, double bestFitness) const override {
        if (currentGeneration < 0) {
            LOG_WARNING("Invalid generation number: " + std::to_string(currentGeneration));
            return false;
        }

        if (!std::isfinite(bestFitness)) {
            LOG_WARNING("Invalid fitness value in shouldStop");
            return false;
        }

        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            if (bestFitness > m_bestFitness) {
                m_bestFitness = bestFitness;
                m_generationsWithoutImprovement = 0;
                LOG_DEBUG("New best fitness found: " + std::to_string(bestFitness));
                return false;
            }
            
            m_generationsWithoutImprovement++;
            LOG_DEBUG("Generations without improvement: " + 
                     std::to_string(m_generationsWithoutImprovement) + 
                     " / " + std::to_string(m_maxGenerationsWithoutImprovement));

            return m_generationsWithoutImprovement >= m_maxGenerationsWithoutImprovement;
        } catch (const std::exception& e) {
            LOG_ERROR("Error in NoImprovementStopping: " + std::string(e.what()));
            return true;  // Stop if there's an error
        }
    }

    void reset() override {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_bestFitness = -std::numeric_limits<double>::infinity();
        m_generationsWithoutImprovement = 0;
        LOG_INFO("NoImprovementStopping criteria reset");
    }

private:
    const int m_maxGenerationsWithoutImprovement;
    mutable double m_bestFitness;
    mutable int m_generationsWithoutImprovement;
    mutable std::mutex m_mutex;
};

class MaxGenerationsStopping : public IStopping {
public:
    explicit MaxGenerationsStopping(int maxGenerations)
        : m_maxGenerations(std::max(1, maxGenerations))
    {
        if (maxGenerations <= 0) {
            LOG_WARNING("Invalid maxGenerations value: " + 
                       std::to_string(maxGenerations) + 
                       ". Using default value of 100");
            m_maxGenerations = 100;
        }
        LOG_INFO("MaxGenerationsStopping initialized with maxGenerations = " + 
                 std::to_string(m_maxGenerations));
    }

    bool shouldStop(int currentGeneration, double bestFitness) const override {
        if (currentGeneration < 0) {
            LOG_WARNING("Invalid generation number: " + std::to_string(currentGeneration));
            return false;
        }

        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            bool shouldStop = currentGeneration >= m_maxGenerations;
            if (shouldStop) {
                LOG_INFO("Maximum number of generations reached: " + 
                        std::to_string(currentGeneration) + 
                        " >= " + std::to_string(m_maxGenerations));
            }
            return shouldStop;
        } catch (const std::exception& e) {
            LOG_ERROR("Error in MaxGenerationsStopping: " + std::string(e.what()));
            return true;  // Stop if there's an error
        }
    }

    void reset() override {
        // Nothing to reset for MaxGenerationsStopping
    }

private:
    const int m_maxGenerations;
    mutable std::mutex m_mutex;
};

class TimeLimitStopping : public IStopping {
public:
    explicit TimeLimitStopping(int limitSeconds)
        : m_limitSeconds(std::max(1, limitSeconds))
    {
        if (limitSeconds <= 0) {
            LOG_WARNING("Invalid time limit: " + 
                       std::to_string(limitSeconds) + 
                       ". Using default value of 3600 seconds (1 hour)");
            m_limitSeconds = 3600;
        }
        reset();  // Initialize m_start
        LOG_INFO("TimeLimitStopping initialized with time limit = " + 
                 std::to_string(m_limitSeconds) + " seconds");
    }

    bool shouldStop(int currentGeneration, double bestFitness) const override {
        if (currentGeneration < 0) {
            LOG_WARNING("Invalid generation number: " + std::to_string(currentGeneration));
            return false;
        }

        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            auto now = std::chrono::steady_clock::now();
            auto elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(
                now - m_start).count();
            
            bool shouldStop = elapsedSeconds >= m_limitSeconds;
            if (shouldStop) {
                LOG_INFO("Time limit reached: " + 
                        std::to_string(elapsedSeconds) + 
                        " >= " + std::to_string(m_limitSeconds) + " seconds");
            }
            return shouldStop;
        } catch (const std::exception& e) {
            LOG_ERROR("Error in TimeLimitStopping: " + std::string(e.what()));
            return true;  // Stop if there's an error
        }
    }

    void reset() override {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_start = std::chrono::steady_clock::now();
        LOG_INFO("TimeLimitStopping criteria reset");
    }

private:
    const int m_limitSeconds;
    mutable std::chrono::steady_clock::time_point m_start;
    mutable std::mutex m_mutex;
}; 