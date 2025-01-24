#include "metaheuristics/path_analyzer.h"
#include "metaheuristics/fitness.h"

PathAnalysisResult PathAnalyzer::analyzePath(
    const std::vector<int>& path,
    const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix
) {
    PathAnalysisResult result;
    std::unordered_set<int> uniqueNodes;
    
    LOG_DEBUG("Starting path analysis...");
    
    for (size_t i = 0; i < path.size() - 1; i++) {
        int current = path[i];
        int next = path[i + 1];
        
        if (current < 0 || current >= (int)adjacencyMatrix.size() || 
            next < 0 || next >= (int)adjacencyMatrix.size()) {
            LOG_ERROR("Invalid vertex index in path");
            continue;
        }
        
        uniqueNodes.insert(current);
        
        const auto& edge = adjacencyMatrix[current][next];
        if (edge.exists) {
            if (edge.weight == 1) {
                result.edgesWeight1++;
            } else {
                result.edgesWeight2or3++;
            }
        }
    }
    
    // Add last node
    if (!path.empty()) {
        uniqueNodes.insert(path.back());
    }
    
    result.uniqueNodesUsed = uniqueNodes.size();
    return result;
} 