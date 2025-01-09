//
// Created by konrad_guest on 07/01/2025.
//



// ================== SimpleFitness ==================
#include <unordered_set>
#include <iostream>

#include "metaheuristics/fitness.h"


double SimpleFitness::evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // decode
    std::string dna = representation->decodeToDNA(individual, instance);
    int k = instance.getK();
    if(k <= 0) return 0.0;

    // build set from dna's k-mers
    std::unordered_set<std::string> kmerSet;
    for(int i=0; i+(k) <= (int)dna.size(); i++){
        kmerSet.insert(dna.substr(i, k));
    }

    // compare with original instance's spectrum
    const auto &spec = instance.getSpectrum();
    std::unordered_set<std::string> spectrumSet(spec.begin(), spec.end());

    // count how many from kmerSet are in spectrum
    int matches = 0;
    for(auto &ss : kmerSet) {
        if(spectrumSet.count(ss)) {
            matches++;
        }
    }

    return (double)matches;
}

double BetterFitness::evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // decode
    std::string dna = representation->decodeToDNA(individual, instance);

    // count matched k-mers
    int k = instance.getK();
    if (k <= 0) return 0.0;

    std::unordered_set<std::string> kmerSet;
    for(int i=0; i + k <= (int)dna.size(); i++){
        kmerSet.insert(dna.substr(i, k));
    }

    const auto &spec = instance.getSpectrum();
    std::unordered_set<std::string> spectrumSet(spec.begin(), spec.end());

    int matches = 0;
    for(auto &ss : kmerSet) {
        if(spectrumSet.count(ss)) {
            matches++;
        }
    }

    // let's incorporate Levenshtein
    std::string originalDNA = instance.getDNA();
    if (originalDNA.empty()) {
        return (double)matches;
    }

    int distance = levenshteinDistance(originalDNA, dna);

    // synergy measure
    double alpha = 0.5; // adjust as needed
    double fitnessVal = matches - alpha * distance;
    return fitnessVal;
}

double SmithWatermanFitness::evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // Get the DNA sequences
    std::string reconstructed = representation->decodeToDNA(individual, instance);
    std::string original = instance.getDNA();
        
    if (reconstructed.empty() || original.empty()) {
        return 0.0;
    }
        
    // Calculate Smith-Waterman score
    int swScore = smithWaterman(original, reconstructed);
        
    // Divide by 5 as mentioned in the paper for comparison purposes
    double normalizedScore = swScore / 5.0;
        
    return normalizedScore;
}

int SmithWatermanFitness::smithWaterman(const std::string& seq1, const std::string& seq2) const {
    int m = seq1.length();
    int n = seq2.length();
    
    // Use single vector instead of 2D
    std::vector<int> score((m+1), 0);
    int maxScore = 0;
    
    for (int j = 1; j <= n; j++) {
        int prev = 0;
        for (int i = 1; i <= m; i++) {
            int temp = score[i];
            int match = prev + (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            int delete_gap = score[i] + GAP_PENALTY;
            int insert_gap = score[i-1] + GAP_PENALTY;
            
            score[i] = std::max({0, match, delete_gap, insert_gap});
            maxScore = std::max(maxScore, score[i]);
            prev = temp;
        }
    }
    return maxScore;
}

std::vector<std::vector<GraphBasedFitness::Edge>> GraphBasedFitness::buildSpectrumGraph(
    const std::vector<std::string>& spectrum, int k) const {
    int n = spectrum.size();
    std::vector<std::vector<Edge>> graph(n);
    
    // Create edges based on overlap between k-mers
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) continue;
            
            int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
            if (weight > 0) {
                graph[i].emplace_back(j, weight);
            }
        }
    }
    
    return graph;
}

int GraphBasedFitness::calculateEdgeWeight(const std::string& from, const std::string& to, int k) const {
    // For k-1 overlap, weight = 1 (perfect overlap)
    // For k-2 overlap, weight = 2 (potential negative error)
    // For k-3 overlap, weight = 3 (likely error)
    
    for (int overlap = k-1; overlap >= k-3; overlap--) {
        if (overlap <= 0) break;
        
        std::string suffix = from.substr(from.length() - overlap);
        std::string prefix = to.substr(0, overlap);
        
        if (suffix == prefix) {
            return k - overlap;  // Weight increases as overlap decreases
        }
    }
    
    return 0;  // No valid overlap found
}

GraphBasedFitness::PathAnalysis GraphBasedFitness::analyzePath(
    const std::vector<int>& path, 
    const std::vector<std::vector<Edge>>& graph) const {
    PathAnalysis analysis;
    
    // Analyze edges in path
    for (size_t i = 0; i < path.size() - 1; i++) {
        int from = path[i];
        int to = path[i + 1];
        
        // Count node usage
        analysis.nodeUsageCount[from]++;
        if (analysis.nodeUsageCount[from] == 1) {
            analysis.uniqueNodesUsed++;
        } else {
            analysis.repeatNodeUsages++;
        }
        
        // Find edge weight
        for (const Edge& edge : graph[from]) {
            if (edge.to == to) {
                if (edge.weight == 1) {
                    analysis.edgesWeight1++;
                } else {
                    analysis.edgesWeight2or3++;
                }
                break;
            }
        }
    }
    
    // Count last node
    analysis.nodeUsageCount[path.back()]++;
    if (analysis.nodeUsageCount[path.back()] == 1) {
        analysis.uniqueNodesUsed++;
    } else {
        analysis.repeatNodeUsages++;
    }
    
    return analysis;
}

std::vector<int> GraphBasedFitness::permutationToPath(void* individual) const {
    auto* perm = static_cast<std::vector<int>*>(individual);
    return *perm;  // Direct conversion as our permutation already represents path
}

double GraphBasedFitness::evaluate(
    void* individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const {
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    if (spectrum.empty() || k <= 0) return 0.0;
    
    // Build spectrum graph
    auto graph = buildSpectrumGraph(spectrum, k);
    
    // Get path from individual
    std::vector<int> path = permutationToPath(individual);
    
    // Analyze path
    PathAnalysis analysis = analyzePath(path, graph);
    
    // Calculate fitness components
    double edgeScore = analysis.edgesWeight1 - (2.0 * analysis.edgesWeight2or3);
    double coverageScore = analysis.uniqueNodesUsed - (0.5 * analysis.repeatNodeUsages);
    
    // Combine scores with weights
    double alpha = 0.7;  // Weight for edge score
    double beta = 0.3;   // Weight for coverage score
    
    return (alpha * edgeScore) + (beta * coverageScore);
}
