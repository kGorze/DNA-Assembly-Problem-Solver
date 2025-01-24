#pragma once

#include "../interfaces/i_crossover.h"
#include "../dna/dna_instance.h"
#include <vector>
#include <memory>

class OnePointCrossover : public ICrossover {
public:
    explicit OnePointCrossover(double crossoverRate) : m_crossoverRate(crossoverRate) {}
    
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

protected:
    std::pair<std::shared_ptr<std::vector<int>>, std::shared_ptr<std::vector<int>>> 
    performCrossover(const std::shared_ptr<std::vector<int>>& parent1,
                    const std::shared_ptr<std::vector<int>>& parent2,
                    const DNAInstance& instance);
private:
    double m_crossoverRate;
};

class OrderCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

protected:
    std::vector<int> performOrderCrossover(
        const std::vector<int>& parent1,
        const std::vector<int>& parent2,
        size_t size);
};

class EdgeRecombination : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
};

class PMXCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
};

class DistancePreservingCrossover : public ICrossover {
private:
    class DistanceMatrix {
    public:
        explicit DistanceMatrix(const std::vector<int>& perm);
        int getDistance(int from, int to) const;
    private:
        std::vector<int> distances;
    };

public:
    std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
}; 