//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef CROSSOVER_H
#define CROSSOVER_H

#include <vector>
#include <memory>
#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"

/**
 * Każdy osobnik to std::shared_ptr<std::vector<int>>.
 */
class ICrossover {
public:
    virtual ~ICrossover() = default;

    virtual std::vector<std::shared_ptr<std::vector<int>>>  
    crossover(const std::vector<std::shared_ptr<std::vector<int>>>& parents,
              const DNAInstance &instance,
              std::shared_ptr<IRepresentation> representation) = 0;
};

// Różne krzyżowania
class OnePointCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> 
    crossover(const std::vector<std::shared_ptr<std::vector<int>>>& parents,
              const DNAInstance &instance,
              std::shared_ptr<IRepresentation> representation) override;
};

// Order Crossover (OX)
class OrderCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> 
    crossover(const std::vector<std::shared_ptr<std::vector<int>>>& parents,
              const DNAInstance &instance,
              std::shared_ptr<IRepresentation> representation) override;
};

// Edge Recombination Crossover (ERX)
class EdgeRecombination : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> 
    crossover(const std::vector<std::shared_ptr<std::vector<int>>>& parents,
              const DNAInstance &instance,
              std::shared_ptr<IRepresentation> representation) override;
};

// Partially Mapped Crossover (PMX)
class PMXCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> 
    crossover(const std::vector<std::shared_ptr<std::vector<int>>>& parents,
              const DNAInstance &instance,
              std::shared_ptr<IRepresentation> representation) override;
};

#endif //CROSSOVER_H
