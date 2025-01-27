#include "../../include/metaheuristics/edge_table.h"
#include "../../include/metaheuristics/dna_utils.h"
#include "../../include/utils/logging.h"
#include <algorithm>

EdgeTable::EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents, const DNAInstance& instance) {
    if (parents.empty()) return;
    
    // Initialize edge lists for each node with overlap quality
    for (const auto& parent : parents) {
        if (!parent) continue;
        
        const auto& genes = parent->getGenes();
        if (genes.empty()) continue;
        
        const auto& spectrum = instance.getSpectrum();
        int k = instance.getK();
        
        for (size_t i = 0; i < genes.size(); ++i) {
            int current = genes[i];
            int prev = i > 0 ? genes[i - 1] : genes.back();
            int next = i < genes.size() - 1 ? genes[i + 1] : genes.front();
            
            // Only add edges if both nodes represent valid k-mers
            if (current >= 0 && static_cast<size_t>(current) < spectrum.size() &&
                prev >= 0 && static_cast<size_t>(prev) < spectrum.size()) {
                // Calculate overlap quality for prev -> current
                int overlapQuality = dna_utils::calculateEdgeWeight(spectrum[prev], spectrum[current], k);
                if (overlapQuality > 0) {
                    edges[current].push_back({prev, overlapQuality});
                }
            }
            
            if (current >= 0 && static_cast<size_t>(current) < spectrum.size() &&
                next >= 0 && static_cast<size_t>(next) < spectrum.size()) {
                // Calculate overlap quality for current -> next
                int overlapQuality = dna_utils::calculateEdgeWeight(spectrum[current], spectrum[next], k);
                if (overlapQuality > 0) {
                    edges[current].push_back({next, overlapQuality});
                }
            }
        }
    }
    
    // Sort edges by overlap quality and remove duplicates
    for (auto& [node, neighbors] : edges) {
        std::sort(neighbors.begin(), neighbors.end(),
                 [](const EdgeInfo& a, const EdgeInfo& b) {
                     return a.overlapQuality > b.overlapQuality;
                 });
        
        // Keep only the highest quality edge for each neighbor
        auto it = std::unique(neighbors.begin(), neighbors.end(),
                            [](const EdgeInfo& a, const EdgeInfo& b) {
                                return a.node == b.node;
                            });
        neighbors.erase(it, neighbors.end());
    }
}

std::vector<EdgeTable::EdgeInfo> EdgeTable::getNeighbors(int node) const {
    auto it = edges.find(node);
    return it != edges.end() ? it->second : std::vector<EdgeInfo>();
}

void EdgeTable::removeNode(int node) {
    // Remove node from all neighbor lists
    for (auto& [_, neighbors] : edges) {
        auto it = std::remove_if(neighbors.begin(), neighbors.end(),
                                 [node](const EdgeInfo& info) {
                                     return info.node == node;
                                 });
        neighbors.erase(it, neighbors.end());
    }
    // Remove node's edge list
    edges.erase(node);
}

bool EdgeTable::hasNode(int node) const {
    return edges.find(node) != edges.end();
} 