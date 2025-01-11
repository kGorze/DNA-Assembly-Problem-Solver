//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef REPLACEMENT_H
#define REPLACEMENT_H

#include <vector>
#include <memory>
#include "metaheuristics/fitness.h"
#include "metaheuristics/population_cache.h"

class IReplacement {
public:
    virtual ~IReplacement() = default;
    virtual std::vector<std::shared_ptr<std::vector<int>>>
    replace(const std::vector<std::shared_ptr<std::vector<int>>> &oldPop,
            const std::vector<std::shared_ptr<std::vector<int>>> &offspring,
            const DNAInstance &instance,
            std::shared_ptr<IFitness> fitness,
            std::shared_ptr<IRepresentation> representation) = 0;
};

class FullReplacement : public IReplacement {
public:
    std::vector<std::shared_ptr<std::vector<int>>>
    replace(const std::vector<std::shared_ptr<std::vector<int>>> &oldPop,
            const std::vector<std::shared_ptr<std::vector<int>>> &offspring,
            const DNAInstance &instance,
            std::shared_ptr<IFitness> fitness,
            std::shared_ptr<IRepresentation> representation) override;
};

class PartialReplacement : public IReplacement {
private:
    double m_replacementRatio;
    std::shared_ptr<IPopulationCache> m_fitnessCache;

public:
    explicit PartialReplacement(double replacementRatio = 0.7, 
                                std::shared_ptr<IPopulationCache> cache = nullptr);
    ~PartialReplacement();
    
    std::vector<std::shared_ptr<std::vector<int>>> replace(
        const std::vector<std::shared_ptr<std::vector<int>>>& oldPop,
        const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) override;
};

#endif //REPLACEMENT_H
