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

Individual::Individual(std::vector<int> genes) : m_genes(std::move(genes)), m_fitness(0.0) {
    validateGenesVector(m_genes);
    m_isValid = true;
}

Individual::Individual(const Individual& other) {
    m_genes = other.m_genes;
    m_fitness = other.m_fitness;
    m_isValid = other.m_isValid;
}

Individual::Individual(Individual&& other) noexcept {
    m_genes = std::move(other.m_genes);
    m_fitness = other.m_fitness;
    m_isValid = other.m_isValid;
    
    // Reset the moved-from object
    other.m_fitness = 0.0;
    other.m_isValid = false;
}

Individual& Individual::operator=(const Individual& other) {
    if (this != &other) {
        m_genes = other.m_genes;
        m_fitness = other.m_fitness;
        m_isValid = other.m_isValid;
    }
    return *this;
}

Individual& Individual::operator=(Individual&& other) noexcept {
    if (this != &other) {
        m_genes = std::move(other.m_genes);
        m_fitness = other.m_fitness;
        m_isValid = other.m_isValid;
        
        // Reset the moved-from object
        other.m_fitness = 0.0;
        other.m_isValid = false;
    }
    return *this;
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

bool Individual::isValid() const {
    return m_isValid && !m_genes.empty();
}

std::string Individual::toString() const {
    std::ostringstream oss;
    oss << "Individual(genes=[";
    for (size_t i = 0; i < m_genes.size(); ++i) {
        if (i > 0) oss << ", ";
        oss << m_genes[i];
    }
    oss << "], fitness=" << m_fitness << ", valid=" << std::boolalpha << m_isValid << ")";
    return oss.str();
}

void Individual::validateGenes() {
    validateGenesVector(m_genes);
    m_isValid = true;
}

void Individual::validateGenesVector(const std::vector<int>& genes) {
    if (genes.empty()) {
        throw std::invalid_argument("Genes vector cannot be empty");
    }
    
    // Only check for negative values, as the upper bound will be checked by the representation
    if (std::any_of(genes.begin(), genes.end(), [](int gene) {
        return gene < 0;
    })) {
        throw std::invalid_argument("Genes vector contains negative values");
    }
}

void Individual::mutate(size_t pos1, size_t pos2) {
    if (pos1 >= m_genes.size() || pos2 >= m_genes.size()) {
        throw std::out_of_range("Mutation positions out of range");
    }
    std::swap(m_genes[pos1], m_genes[pos2]);
}

void Individual::reverse(size_t start, size_t end) {
    if (start >= m_genes.size() || end >= m_genes.size() || start > end) {
        throw std::out_of_range("Invalid reverse range");
    }
    std::reverse(m_genes.begin() + start, m_genes.begin() + end + 1);
}

void Individual::shift(size_t start, size_t end, int positions) {
    if (start >= m_genes.size() || end >= m_genes.size() || start > end) {
        throw std::out_of_range("Invalid shift range");
    }
    
    if (positions == 0) return;
    
    std::vector<int> temp(m_genes.begin() + start, m_genes.begin() + end + 1);
    size_t length = end - start + 1;
    positions = ((positions % length) + length) % length; // Normalize positions
    
    std::rotate(temp.begin(), temp.begin() + positions, temp.end());
    std::copy(temp.begin(), temp.end(), m_genes.begin() + start);
} 