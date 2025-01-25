#include "../../include/metaheuristics/individual.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <memory>

// Currently all methods are implemented inline in the header file
// This file exists for future implementations of more complex methods 

Individual::Individual(std::vector<int> genes) 
    : m_genes(std::move(genes)) 
{
    validateGenes();
}

Individual::Individual(const Individual& other) 
    : m_genes(other.m_genes)
    , m_fitness(other.m_fitness)
    , m_isValid(other.m_isValid) 
{
    validateGenes();  // Revalidate after copy
}

Individual::Individual(Individual&& other) noexcept 
    : m_genes(std::move(other.m_genes))
    , m_fitness(other.m_fitness)
    , m_isValid(other.m_isValid)
{
    // Clear source object
    other.m_genes.clear();
    other.m_fitness = 0.0;
    other.m_isValid = false;
}

Individual& Individual::operator=(const Individual& other) {
    if (this != &other) {
        try {
            std::vector<int> temp(other.m_genes);  // Copy into temporary
            m_genes = std::move(temp);             // Move into member
            m_fitness = other.m_fitness;
            m_isValid = other.m_isValid;
            validateGenes();  // Revalidate after assignment
        } catch (const std::exception& e) {
            // If assignment fails, ensure object is in valid state
            m_genes.clear();
            m_fitness = 0.0;
            m_isValid = false;
            throw;  // Rethrow the exception
        }
    }
    return *this;
}

Individual& Individual::operator=(Individual&& other) noexcept {
    if (this != &other) {
        m_genes = std::move(other.m_genes);
        m_fitness = other.m_fitness;
        m_isValid = other.m_isValid;
        
        // Clear source object
        other.m_genes.clear();
        other.m_fitness = 0.0;
        other.m_isValid = false;
    }
    return *this;
}

void Individual::setGenes(std::vector<int> genes) {
    try {
        std::vector<int> temp = std::move(genes);  // Move into temporary
        validateGenesVector(temp);                 // Validate before assigning
        m_genes = std::move(temp);                 // Move into member
        m_isValid = true;
    } catch (const std::exception& e) {
        m_isValid = false;
        throw;  // Rethrow the exception
    }
}

void Individual::setFitness(double fitness) {
    if (!std::isfinite(fitness)) {
        throw std::invalid_argument("Fitness must be a finite number");
    }
    if (fitness < -std::numeric_limits<double>::max() || 
        fitness > std::numeric_limits<double>::max()) {
        throw std::invalid_argument("Fitness value out of range");
    }
    m_fitness = fitness;
}

std::string Individual::toString() const {
    if (empty()) {
        return "empty";
    }
    
    try {
        std::stringstream ss;
        ss.exceptions(std::ios::badbit | std::ios::failbit);
        
        for (int gene : m_genes) {
            ss << gene << " ";
        }
        
        std::string result = ss.str();
        if (!result.empty()) {
            result.pop_back();  // Remove trailing space
        }
        return result;
    } catch (const std::exception& e) {
        return "error: " + std::string(e.what());
    }
}

void Individual::validateGenes() {
    try {
        validateGenesVector(m_genes);
    } catch (const std::exception&) {
        m_isValid = false;
        throw;
    }
}

void Individual::validateGenesVector(const std::vector<int>& genes) {
    if (genes.empty()) {
        throw std::invalid_argument("Gene sequence cannot be empty");
    }
    
    // Check for invalid values
    auto invalidIt = std::find_if(genes.begin(), genes.end(),
        [](int gene) { 
            return gene < 0 || gene > std::numeric_limits<int>::max() / 2; 
        });
    
    if (invalidIt != genes.end()) {
        throw std::invalid_argument(
            "Invalid gene value: " + std::to_string(*invalidIt));
    }
    
    // Check for duplicate genes if needed
    auto sorted = genes;
    std::sort(sorted.begin(), sorted.end());
    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end()) {
        throw std::invalid_argument("Duplicate genes detected");
    }
    
    m_isValid = true;
} 