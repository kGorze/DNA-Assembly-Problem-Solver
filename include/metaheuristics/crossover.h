//
// Created by konrad_guest on 07/01/2025.
//

#ifndef CROSSOVER_H
#define CROSSOVER_H

#include <vector>
#include "metaheuristics/crossover.h"
#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"

class ICrossover {
public:
    virtual ~ICrossover() = default;
    /**
     * Metoda przyjmuje rodziców i zwraca potomstwo
     */
    //virtual std::vector<std::vector<double>> crossover(const std::vector<std::vector<double>> &parents) = 0;

    virtual std::vector<void*>  crossover(const std::vector<void*> &parents,const DNAInstance &instance, std::shared_ptr<IRepresentation> representation) = 0;
};

// Różne krzyżowania
class OnePointCrossover : public ICrossover {
public:
    std::vector<void*>  crossover(const std::vector<void*> &parents,const DNAInstance &instance, std::shared_ptr<IRepresentation> representation) override;
};

// class TwoPointCrossover : public ICrossover {
// public:
//     std::vector<std::vector<double>> 
//     crossover(const std::vector<std::vector<double>> &parents) override;
// };
//
// class UniformCrossover : public ICrossover {
// public:
//     std::vector<std::vector<double>> 
//     crossover(const std::vector<std::vector<double>> &parents) override;
// };
//
// class ArithmeticCrossover : public ICrossover {
// public:
//     std::vector<std::vector<double>> 
//     crossover(const std::vector<std::vector<double>> &parents) override;
// };


#endif //CROSSOVER_H
