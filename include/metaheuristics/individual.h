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
    Individual() = default;
    
    // Delete copy operations due to mutex member
    Individual(const Individual&) = delete;
    Individual& operator=(const Individual&) = delete;
    
    // Allow move operations
    Individual(Individual&&) noexcept = default;
    Individual& operator=(Individual&&) noexcept = default;
    
    virtual ~Individual() = default;

    // Getters and setters with validation
    const std::vector<int>& getGenes() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_genes;
    }
    
    std::vector<int>& getGenes() {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_genes;
    }
    
    void setGenes(const std::vector<int>& genes) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_genes = genes;
    }
    
    void setGenes(std::vector<int>&& genes) {
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
    double m_fitness = std::numeric_limits<double>::infinity();
    bool m_isValid = false;
    mutable std::mutex m_mutex;  // For thread safety
}; 