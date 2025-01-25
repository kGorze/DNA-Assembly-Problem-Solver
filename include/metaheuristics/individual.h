#pragma once

#include <vector>
#include <string>
#include <memory>
#include <sstream>

class Individual {
public:
    // Constructors
    Individual() = default;
    
    explicit Individual(std::vector<int> genes) : m_genes(std::move(genes)) {}
    
    Individual(const Individual& other) : m_genes(other.m_genes), m_fitness(other.m_fitness) {}
    
    Individual(Individual&& other) noexcept 
        : m_genes(std::move(other.m_genes)), m_fitness(other.m_fitness) {}

    // Assignment operators
    Individual& operator=(const Individual& other) {
        if (this != &other) {
            m_genes = other.m_genes;
            m_fitness = other.m_fitness;
        }
        return *this;
    }
    
    Individual& operator=(Individual&& other) noexcept {
        if (this != &other) {
            m_genes = std::move(other.m_genes);
            m_fitness = other.m_fitness;
        }
        return *this;
    }

    // Getters and setters
    const std::vector<int>& getGenes() const { return m_genes; }
    std::vector<int>& getGenes() { return m_genes; }
    void setGenes(std::vector<int> genes) { m_genes = std::move(genes); }
    
    double getFitness() const { return m_fitness; }
    void setFitness(double fitness) { m_fitness = fitness; }

    // Utility methods
    bool empty() const { return m_genes.empty(); }
    size_t size() const { return m_genes.size(); }
    
    std::string toString() const {
        std::stringstream ss;
        for (int gene : m_genes) {
            ss << gene;
        }
        return ss.str();
    }

private:
    std::vector<int> m_genes;
    double m_fitness = 0.0;
}; 