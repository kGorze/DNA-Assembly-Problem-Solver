#pragma once

#include <vector>
#include <string>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <mutex>
#include <cmath>

/**
 * Class representing an individual in the genetic algorithm.
 */
class Individual {
public:
    Individual() = default;
    explicit Individual(std::vector<int> genes);
    Individual(const Individual& other);
    Individual(Individual&& other) noexcept;
    ~Individual() = default;

    Individual& operator=(const Individual& other);
    Individual& operator=(Individual&& other) noexcept;

    // Getters
    const std::vector<int>& getGenes() const { return m_genes; }
    std::vector<int>& getGenes() { return m_genes; }  // Non-const version for mutation
    double getFitness() const { return m_fitness; }
    size_t getSize() const { return m_genes.size(); }

    // Setters
    void setGenes(const std::vector<int>& genes) {
        m_genes = genes;
        m_isValid = !m_genes.empty();  // Basic validation - genes vector is not empty
    }
    void setFitness(double fitness);

    // Utility methods
    bool isValid() const;
    std::string toString() const;
    void mutate(size_t pos1, size_t pos2);
    void reverse(size_t start, size_t end);
    void shift(size_t start, size_t end, int positions);

    // Validation methods
    void validate() {
        if (!m_isValid) {
            validateGenes();
        }
    }

    std::shared_ptr<Individual> clone() const {
        auto newIndividual = std::make_shared<Individual>();
        newIndividual->setGenes(m_genes);
        newIndividual->setFitness(m_fitness);
        return newIndividual;
    }

private:
    void validateGenes();
    void validateGenesVector(const std::vector<int>& genes);

    std::vector<int> m_genes;
    double m_fitness = 0.0;
    bool m_isValid = false;
    mutable std::mutex m_mutex;  // For thread safety
}; 