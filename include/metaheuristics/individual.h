#pragma once

#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <mutex>

class Individual {
public:
    // Constructors
    Individual() = default;
    explicit Individual(std::vector<int> genes);
    Individual(const Individual& other);
    Individual(Individual&& other) noexcept;

    // Assignment operators
    Individual& operator=(const Individual& other);
    Individual& operator=(Individual&& other) noexcept;

    // Destructor
    ~Individual() = default;

    // Getters and setters with validation
    const std::vector<int>& getGenes() const { 
        if (!m_isValid) {
            throw std::runtime_error("Individual is in an invalid state");
        }
        return m_genes; 
    }
    
    std::vector<int>& getGenes() { 
        m_isValid = false;  // Mark as needing validation since genes can be modified
        return m_genes; 
    }
    
    void setGenes(std::vector<int> genes);
    
    double getFitness() const { 
        if (!m_isValid) {
            throw std::runtime_error("Individual is in an invalid state");
        }
        return m_fitness; 
    }
    
    void setFitness(double fitness) {
        if (!m_isValid) {
            throw std::runtime_error("Individual is in an invalid state");
        }
        if (!std::isfinite(fitness)) {
            throw std::invalid_argument("Fitness must be a finite number");
        }
        m_fitness = fitness;
    }

    // Utility methods
    bool empty() const { return m_genes.empty(); }
    size_t size() const { return m_genes.size(); }
    bool isValid() const { return m_isValid; }
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
    double m_fitness{0.0};
    bool m_isValid{false};
    mutable std::mutex m_mutex;  // For thread safety
}; 