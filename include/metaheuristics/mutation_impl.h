#pragma once

#include "../interfaces/i_mutation.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include <random>
#include <algorithm>

class PointMutation : public IMutation {
public:
    explicit PointMutation(double mutationRate) 
        : m_mutationRate(mutationRate)
    {
        if (m_mutationRate < 0.0 || m_mutationRate > 1.0) {
            LOG_WARNING("Invalid mutation rate " + std::to_string(m_mutationRate) + " - clamping to valid range");
            m_mutationRate = std::clamp(m_mutationRate, 0.0, 1.0);
        }
    }

    void mutate(
        std::shared_ptr<std::vector<int>>& solution,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

    double getMutationRate() const { return m_mutationRate; }

private:
    double m_mutationRate;
}; 