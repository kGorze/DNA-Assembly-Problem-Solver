//
// Created by konrad_guest on 07/01/2025.
// SMART

#pragma once

#include <vector>
#include <memory>
#include "generator/dna_generator.h"

// Forward declarations
class IFitness;
class IRepresentation;
class GAConfig;  // Forward declare GAConfig instead of including it

class IStopping {
public:
    virtual ~IStopping() = default;
    virtual bool stop(const std::vector<std::shared_ptr<std::vector<int>>>& population,
                     int generation,
                     const DNAInstance& instance,
                     std::shared_ptr<IFitness> fitness,
                     std::shared_ptr<IRepresentation> representation) = 0;
};

class MaxGenerationsStopping : public IStopping {
public:
    // Constructor taking GAConfig reference
    explicit MaxGenerationsStopping(GAConfig& config);
    
    // Constructor with explicit max generations value
    explicit MaxGenerationsStopping(int maxGen);
    
    bool stop(const std::vector<std::shared_ptr<std::vector<int>>>& population,
              int generation,
              const DNAInstance& instance,
              std::shared_ptr<IFitness> fitness,
              std::shared_ptr<IRepresentation> representation) override;

    // Add method to set maxGenerations
    void setMaxGenerations(int maxGen) {
        m_maxGenerations = maxGen;
        m_useConfig = false;
        std::cout << "[MaxGenerationsStopping] Set fixed maxGenerations = " << maxGen << std::endl;
    }

private:
    int m_maxGenerations;  // Store locally if provided in constructor
    bool m_useConfig;      // Whether to use GAConfig or local value
};

class NoImprovementStopping : public IStopping {
private:
    double m_bestFitness;
    int m_generationsWithoutImprovement;
    int m_maxGenerationsWithoutImprovement;
    int m_totalGenerations;
    
public:
    NoImprovementStopping(int maxGenerations) 
        : m_bestFitness(-std::numeric_limits<double>::infinity())
        , m_generationsWithoutImprovement(0)
        , m_totalGenerations(maxGenerations)
        , m_maxGenerationsWithoutImprovement(static_cast<int>(maxGenerations * 0.3)) // 30%
    {}
    bool stop(const std::vector<std::shared_ptr<std::vector<int>>> &population,
              int generation,
              const DNAInstance &instance,
              std::shared_ptr<IFitness> fitness,
              std::shared_ptr<IRepresentation> representation) override;
};
