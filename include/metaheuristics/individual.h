#pragma once

#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <mutex>
#include <cmath>

class Individual {
public:
    // Constructors
    Individual() : m_fitness(0.0) {}
    explicit Individual(std::vector<int> genes) : m_genes(std::move(genes)), m_fitness(0.0) {}
    
    // Delete copy operations due to mutex member
    Individual(const Individual&) = delete;
    Individual& operator=(const Individual&) = delete;
    
    // Allow move operations
    Individual(Individual&& other) noexcept {
        std::lock_guard<std::mutex> lock(other.m_mutex);
        m_genes = std::move(other.m_genes);
        m_fitness = other.m_fitness;
    }
    
    Individual& operator=(Individual&& other) noexcept {
        if (this != &other) {
            std::lock_guard<std::mutex> lock1(m_mutex);
            std::lock_guard<std::mutex> lock2(other.m_mutex);
            m_genes = std::move(other.m_genes);
            m_fitness = other.m_fitness;
        }
        return *this;
    }
    
    virtual ~Individual() = default;

    // Getters and setters with validation
    const std::vector<int>& getGenes() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_genes;
    }
    
    void setGenes(std::vector<int> genes) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_genes = std::move(genes);
    }
    
    double getFitness() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_fitness;
    }
    
    void setFitness(double fitness) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (!std::isfinite(fitness)) {
            throw std::invalid_argument("Fitness must be a finite number");
        }
        m_fitness = fitness;
        m_isValid = true;
    }

    // Utility methods
    bool empty() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_genes.empty();
    }
    
    size_t size() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_genes.size();
    }
    
    bool isValid() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_isValid;
    }
    
    void setValid(bool valid) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_isValid = valid;
    }
    
    std::string toString() const;

    // Validation methods
    void validate() {
        if (!m_isValid) {
            validateGenes();
        }
    }

private:
    void validateGenes();
    void validateGenesVector(const std::vector<int>& genes);

    std::vector<int> m_genes;
    double m_fitness;
    bool m_isValid = false;
    mutable std::mutex m_mutex;  // For thread safety
}; 