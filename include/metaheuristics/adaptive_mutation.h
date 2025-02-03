#pragma once

#include "../configuration/genetic_algorithm_configuration.h"
#include <memory>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <mutex>
#include <atomic>

// Forward declaration
struct DummyAdaptiveParams;

class AdaptiveMutation {
private:
    // Thread-safe statistics management
    class OperatorStats {
    private:
        mutable std::mutex m_mutex;
        const size_t m_size;  // Move m_size before vectors
        std::vector<double> m_probabilities;
        std::vector<int> m_usage;

    public:
        explicit OperatorStats(size_t size) 
            : m_size(size)  // Initialize m_size first
            , m_probabilities(size, 1.0/static_cast<double>(size))
            , m_usage(size, 0) {}

        void reset() {
            std::lock_guard<std::mutex> lock(m_mutex);
            std::fill(m_usage.begin(), m_usage.end(), 0);
            std::fill(m_probabilities.begin(), m_probabilities.end(), 
                     1.0/static_cast<double>(m_size));
        }

        void updateStats(size_t index, bool success) {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (index >= m_size) {
                throw std::out_of_range("Invalid operator index");
            }
            
            m_usage[index]++;
            if (success) {
                m_probabilities[index] *= 1.1;
            } else {
                m_probabilities[index] *= 0.9;
            }
            normalize();
        }

        std::vector<double> getProbabilities() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_probabilities;
        }

        std::vector<int> getUsage() const {
            std::lock_guard<std::mutex> lock(m_mutex);
            return m_usage;
        }

    private:
        void normalize() {
            // Called with mutex already locked
            double sum = std::accumulate(m_probabilities.begin(), m_probabilities.end(), 0.0);
            if (sum > 0.0) {
                for (auto& prob : m_probabilities) {
                    prob /= sum;
                }
            } else {
                std::fill(m_probabilities.begin(), m_probabilities.end(), 
                         1.0/static_cast<double>(m_size));
            }
        }
    };

    // Configuration state - reorder members to match initialization order
    const std::shared_ptr<const GAConfig> m_config;
    std::atomic<double> m_currentMutationRate;
    std::atomic<int> m_stagnationCounter;
    std::atomic<double> m_lastBestFitness;
    OperatorStats m_stats;
    
    static constexpr int ADAPTATION_INTERVAL = 100;
    static constexpr size_t NUM_OPERATORS = 3;
    
    mutable std::mutex m_configMutex;

    static void validateConfig(const std::shared_ptr<const GAConfig>& config) {
        if (!config) {
            throw std::invalid_argument("Config cannot be null");
        }
        const auto& params = config->getAdaptiveParams();
        if (params.minMutationRate < 0.0 || params.maxMutationRate > 1.0 || 
            params.minMutationRate >= params.maxMutationRate) {
            throw std::invalid_argument("Invalid mutation rate bounds");
        }
    }

public:
    explicit AdaptiveMutation(std::shared_ptr<const GAConfig> config)
        : m_config(std::move(config))
        , m_currentMutationRate(0.0)
        , m_stagnationCounter(0)
        , m_lastBestFitness(0.0)
        , m_stats(NUM_OPERATORS) {
        validateConfig(m_config);
        m_currentMutationRate = m_config->getMutationRate();
    }
    
    // Prevent copying and moving
    AdaptiveMutation(const AdaptiveMutation&) = delete;
    AdaptiveMutation& operator=(const AdaptiveMutation&) = delete;
    AdaptiveMutation(AdaptiveMutation&&) = delete;
    AdaptiveMutation& operator=(AdaptiveMutation&&) = delete;
    
    void updateMutationRate(double bestFitness) {
        std::lock_guard<std::mutex> lock(m_configMutex);
        
        if (!m_config) {
            throw std::runtime_error("Config is null");
        }
        
        const auto& params = m_config->getAdaptiveParams();
        if (!params.useAdaptiveMutation) return;
        
        double lastFitness = m_lastBestFitness.load();
        double improvement = bestFitness - lastFitness;
        
        if (improvement > params.improvementThreshold) {
            m_stagnationCounter.store(0);
            m_currentMutationRate.store(std::max(
                params.minMutationRate,
                m_currentMutationRate.load() * 0.95
            ));
        } else {
            int stagnation = m_stagnationCounter.fetch_add(1) + 1;
            if (stagnation >= params.stagnationGenerations) {
                m_currentMutationRate.store(std::min(
                    params.maxMutationRate,
                    m_currentMutationRate.load() * 1.1
                ));
                m_stagnationCounter.store(0);
            }
        }
        
        m_lastBestFitness.store(bestFitness);
    }
    
    void updateOperatorStats(size_t operatorIndex, bool success) {
        m_stats.updateStats(operatorIndex, success);
    }
    
    void reset() {
        std::lock_guard<std::mutex> lock(m_configMutex);
        
        if (!m_config) {
            throw std::runtime_error("Config is null");
        }
        
        m_currentMutationRate.store(m_config->getMutationRate());
        m_stagnationCounter.store(0);
        m_lastBestFitness.store(0.0);
        m_stats.reset();
    }
    
    [[nodiscard]] double getCurrentMutationRate() const {
        return m_currentMutationRate.load();
    }
    
    [[nodiscard]] std::vector<double> getOperatorProbabilities() const {
        return m_stats.getProbabilities();
    }
    
    [[nodiscard]] std::vector<int> getOperatorUsage() const {
        return m_stats.getUsage();
    }
    
    // Method to clear any references before destruction
    void clearReferences() {
        std::lock_guard<std::mutex> lock(m_configMutex);
        m_currentMutationRate.store(0.0);
        m_stagnationCounter.store(0);
        m_lastBestFitness.store(0.0);
        m_stats.reset();
        // We don't reset m_config here as it's const and will be handled by the destructor
    }
    
    // Add explicit destructor to ensure proper cleanup order
    ~AdaptiveMutation() = default;
}; 