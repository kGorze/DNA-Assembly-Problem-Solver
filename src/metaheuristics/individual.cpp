#include "../../include/metaheuristics/individual.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <memory>
#include <mutex>
#include <sstream>

// Currently all methods are implemented inline in the header file
// This file exists for future implementations of more complex methods 

Individual::Individual(std::vector<int> genes) : m_genes(std::move(genes)) {
    validateGenesVector(m_genes);
    m_isValid = true;
}

Individual::Individual(const Individual& other) {
    std::lock_guard<std::mutex> lock(other.m_mutex);
    m_genes = other.m_genes;
    m_fitness = other.m_fitness;
    m_isValid = other.m_isValid;
}

Individual::Individual(Individual&& other) noexcept {
    std::lock_guard<std::mutex> lock(other.m_mutex);
    m_genes = std::move(other.m_genes);
    m_fitness = other.m_fitness;
    m_isValid = other.m_isValid;
    
    // Reset the moved-from object
    other.m_fitness = 0.0;
    other.m_isValid = false;
}

Individual& Individual::operator=(const Individual& other) {
    if (this != &other) {
        std::lock_guard<std::mutex> lock1(m_mutex);
        std::lock_guard<std::mutex> lock2(other.m_mutex);
        
        std::vector<int> tempGenes = other.m_genes;  // Copy genes first
        validateGenesVector(tempGenes);  // Validate before modifying state
        
        m_genes = std::move(tempGenes);
        m_fitness = other.m_fitness;
        m_isValid = other.m_isValid;
    }
    return *this;
}

Individual& Individual::operator=(Individual&& other) noexcept {
    if (this != &other) {
        std::lock_guard<std::mutex> lock1(m_mutex);
        std::lock_guard<std::mutex> lock2(other.m_mutex);
        
        m_genes = std::move(other.m_genes);
        m_fitness = other.m_fitness;
        m_isValid = other.m_isValid;
        
        // Reset the moved-from object
        other.m_fitness = 0.0;
        other.m_isValid = false;
    }
    return *this;
}

void Individual::setGenes(std::vector<int> genes) {
    std::lock_guard<std::mutex> lock(m_mutex);
    validateGenesVector(genes);  // Validate before modifying state
    m_genes = std::move(genes);
    m_isValid = true;
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
    std::lock_guard<std::mutex> lock(m_mutex);
    std::ostringstream oss;
    oss << "Individual(genes=[";
    if (!m_genes.empty()) {
        for (size_t i = 0; i < m_genes.size() - 1; ++i) {
            oss << m_genes[i] << ", ";
        }
        oss << m_genes.back();
    }
    oss << "], fitness=" << m_fitness << ", valid=" << std::boolalpha << m_isValid << ")";
    return oss.str();
}

void Individual::validateGenes() {
    std::lock_guard<std::mutex> lock(m_mutex);
    validateGenesVector(m_genes);
    m_isValid = true;
}

void Individual::validateGenesVector(const std::vector<int>& genes) {
    if (genes.empty()) {
        throw std::invalid_argument("Genes vector cannot be empty");
    }
    
    // Check for duplicate genes
    std::vector<int> sortedGenes = genes;
    std::sort(sortedGenes.begin(), sortedGenes.end());
    if (std::adjacent_find(sortedGenes.begin(), sortedGenes.end()) != sortedGenes.end()) {
        throw std::invalid_argument("Genes vector contains duplicate values");
    }
    
    // Check for invalid gene values
    if (std::any_of(genes.begin(), genes.end(), [](int gene) {
        return gene < 0 || gene >= std::numeric_limits<int>::max();
    })) {
        throw std::invalid_argument("Genes vector contains invalid values");
    }
} 