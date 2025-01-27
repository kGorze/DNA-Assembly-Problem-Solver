#pragma once

#include <vector>
#include <memory>
#include <unordered_map>
#include "individual.h"
#include "../dna/dna_instance.h"

class EdgeTable {
public:
    struct EdgeInfo {
        int node;
        int overlapQuality;
        bool operator==(const EdgeInfo& other) const { return node == other.node; }
    };
    
    EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents, const DNAInstance& instance);
    std::vector<EdgeInfo> getNeighbors(int node) const;
    void removeNode(int node);
    bool hasNode(int node) const;
    
private:
    std::unordered_map<int, std::vector<EdgeInfo>> edges;
}; 