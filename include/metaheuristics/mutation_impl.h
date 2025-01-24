#pragma once

#include "../interfaces/i_mutation.h"
#include "../dna/dna_instance.h"

class PointMutation : public IMutation {
public:
    explicit PointMutation(double mutationRate) : m_mutationRate(mutationRate) {}
    
    void mutate(
        std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

private:
    double m_mutationRate;
}; 