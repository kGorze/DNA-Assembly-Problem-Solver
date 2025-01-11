//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef STOPPING_CRITERIA_H
#define STOPPING_CRITERIA_H

#include <vector>
#include <memory>
#include "metaheuristics/fitness.h"
#include "generator/dna_generator.h"
#include "metaheuristics/representation.h"

class IStopping {
public:
    virtual ~IStopping() = default;

    virtual bool stop(const std::vector<std::shared_ptr<std::vector<int>>> &population,
                      int generation,
                      const DNAInstance &instance,
                      std::shared_ptr<IFitness> fitness,
                      std::shared_ptr<IRepresentation> representation) = 0;
};

class MaxGenerationsStopping : public IStopping {
public:
    explicit MaxGenerationsStopping(int maxGen) : m_maxGen(maxGen) {}

    bool stop(const std::vector<std::shared_ptr<std::vector<int>>> &population,
              int generation,
              const DNAInstance &instance,
              std::shared_ptr<IFitness> fitness,
              std::shared_ptr<IRepresentation> representation) override
    {
        return generation >= m_maxGen;
    }

private:
    int m_maxGen;
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

#endif //STOPPING_CRITERIA_H
