#pragma once

#include "../configuration/genetic_algorithm_configuration.h"
#include <memory>
#include <vector>

class AdaptiveMutation {
private:
    const GAConfig& m_config;
    double m_currentMutationRate;
    int m_stagnationCounter;
    double m_lastBestFitness;
    
public:
    explicit AdaptiveMutation(const GAConfig& config) 
        : m_config(config)
        , m_currentMutationRate(config.getMutationRate())
        , m_stagnationCounter(0)
        , m_lastBestFitness(0.0) {}
    
    // Update mutation rate based on population performance
    void updateMutationRate(double bestFitness) {
        const auto& params = m_config.getAdaptiveParams();
        if (!params.useAdaptiveMutation) return;
        
        // Check for improvement
        double improvement = bestFitness - m_lastBestFitness;
        if (improvement > params.improvementThreshold) {
            // Reset stagnation counter on significant improvement
            m_stagnationCounter = 0;
            // Slightly decrease mutation rate since we're improving
            m_currentMutationRate = std::max(
                params.minMutationRate,
                m_currentMutationRate * 0.95
            );
        } else {
            // Increment stagnation counter
            m_stagnationCounter++;
            
            // If stagnated for too long, increase mutation rate
            if (m_stagnationCounter >= params.stagnationGenerations) {
                m_currentMutationRate = std::min(
                    params.maxMutationRate,
                    m_currentMutationRate * 1.1
                );
                m_stagnationCounter = 0;  // Reset counter
            }
        }
        
        m_lastBestFitness = bestFitness;
    }
    
    // Get current mutation rate
    double getCurrentMutationRate() const {
        return m_currentMutationRate;
    }
    
    // Reset adaptation
    void reset() {
        m_currentMutationRate = m_config.getMutationRate();
        m_stagnationCounter = 0;
        m_lastBestFitness = 0.0;
    }
}; 