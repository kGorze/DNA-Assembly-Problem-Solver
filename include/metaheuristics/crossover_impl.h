#pragma once

#include "../interfaces/i_crossover.h"
#include "../dna/dna_instance.h"
#include "individual.h"
#include <vector>
#include <memory>
#include <random>

class OnePointCrossover : public ICrossover {
public:
    explicit OnePointCrossover(double crossoverRate) : m_crossoverRate(crossoverRate) {}
    
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

protected:
    std::pair<std::shared_ptr<Individual>, std::shared_ptr<Individual>> 
    performCrossover(const std::shared_ptr<Individual>& parent1,
                    const std::shared_ptr<Individual>& parent2,
                    const DNAInstance& instance,
                    std::shared_ptr<IRepresentation> representation);
private:
    double m_crossoverRate;
    std::random_device rd;
    std::mt19937 gen{rd()};
};

class OrderCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

protected:
    std::vector<int> performOrderCrossover(
        const std::vector<int>& parent1,
        const std::vector<int>& parent2,
        size_t size);
private:
    std::random_device rd;
    std::mt19937 gen{rd()};
};

class EdgeRecombination : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;
    
protected:
    class EdgeTable {
    public:
        explicit EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents);
        std::vector<int> getNeighbors(int node) const;
        void removeNode(int node);
        bool hasNode(int node) const;
        
    private:
        std::unordered_map<int, std::vector<int>> edges;
    };
    
private:
    std::random_device rd;
    std::mt19937 gen{rd()};
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