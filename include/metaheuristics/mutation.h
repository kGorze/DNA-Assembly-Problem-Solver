//
// Created by konrad_guest on 07/01/2025.
//

#ifndef MUTATION_H
#define MUTATION_H

#include <vector>
#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"

class IMutation {
public:
    virtual ~IMutation() = default;
    virtual void mutate(std::vector<void*>        &offspring,
                           const DNAInstance         &instance,
                           std::shared_ptr<IRepresentation> representation) = 0;
};

// Różne rodzaje mutacji:
class PointMutation : public IMutation {
public:
    PointMutation(double rate = 0.1) : mutationRate(rate) {}
    void setMutationRate(double rate) { mutationRate = rate; }
    double getMutationRate() const { return mutationRate; }
    void mutate(std::vector<void*> &offspring,
                const DNAInstance &instance,
                std::shared_ptr<IRepresentation> representation) override;
private:
    double mutationRate;
};



// class GaussianMutation : public IMutation {
// public:
//     void mutate(std::vector<std::vector<double>> &offspring) override;
// };
//
// class InversionMutation : public IMutation {
// public:
//     void mutate(std::vector<std::vector<double>> &offspring) override;
// };

#endif //MUTATION_H
