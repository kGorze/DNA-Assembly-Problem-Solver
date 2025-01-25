#pragma once

#include <vector>
#include <string>
#include <memory>
#include <sstream>

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
    const std::vector<int>& getGenes() const { return m_genes; }
    
    std::vector<int>& getGenes() { 
        m_isValid = false;  // Mark as needing validation since genes can be modified
        return m_genes; 
    }
    
    void setGenes(std::vector<int> genes);
    double getFitness() const { return m_fitness; }
    void setFitness(double fitness);

    // Utility methods
    bool empty() const { return m_genes.empty(); }
    size_t size() const { return m_genes.size(); }
    bool isValid() const { return m_isValid; }
    std::string toString() const;

    // Validation methods
    void validate() { validateGenes(); }

private:
    void validateGenes();
    void validateGenesVector(const std::vector<int>& genes);

    std::vector<int> m_genes;
    double m_fitness{0.0};
    bool m_isValid{false};
}; 