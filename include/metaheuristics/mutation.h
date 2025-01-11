//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef MUTATION_H
#define MUTATION_H

#include <vector>
#include <memory>
#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"

class IMutation {
public:
    virtual ~IMutation() = default;
    virtual void mutate(std::vector<std::shared_ptr<std::vector<int>>> &offspring,
                        const DNAInstance &instance,
                        std::shared_ptr<IRepresentation> representation) = 0;
};

class PointMutation : public IMutation {
public:
    PointMutation(double rate = 0.1) : mutationRate(rate) {}
    void setMutationRate(double rate) { mutationRate = rate; }
    double getMutationRate() const { return mutationRate; }
    void mutate(std::vector<std::shared_ptr<std::vector<int>>> &offspring,
                const DNAInstance &instance,
                std::shared_ptr<IRepresentation> representation) override;
private:
    double mutationRate;
};

#endif //MUTATION_H
