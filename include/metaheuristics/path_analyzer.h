#ifndef PATH_ANALYZER_H
#define PATH_ANALYZER_H

#include <vector>
#include <unordered_set>
#include "utils/logging.h"
#include "metaheuristics/preprocessed_edge.h"

// Forward declarations
class OptimizedGraphBasedFitness;

struct PathAnalysisResult {
    int uniqueNodesUsed;
    double edgesWeight1;
    double edgesWeight2or3;
    
    PathAnalysisResult() : uniqueNodesUsed(0), edgesWeight1(0.0), edgesWeight2or3(0.0) {}
};

class PathAnalyzer {
public:
    static PathAnalysisResult analyzePath(
        const std::vector<int>& path,
        const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix
    );
};

#endif // PATH_ANALYZER_H 