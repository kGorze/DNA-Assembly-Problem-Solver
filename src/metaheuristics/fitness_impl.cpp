#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <numeric>
#include <cmath>

double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const {
    if (!individual || !representation) {
        return 0.0;
    }

    const auto& genes = individual->getGenes();
    
    // Get DNA string representation for Levenshtein calculation
    std::string dnaStr = representation->toString(individual, instance);
    std::string targetSequence = instance.getTargetSequence();
    
    // Always calculate and log Levenshtein distance if target sequence is available
    if (!targetSequence.empty()) {
        int levenshteinDist = calculateLevenshteinDistance(dnaStr, targetSequence);
        double levenshteinRatio = static_cast<double>(levenshteinDist) / targetSequence.length();
        LOG_INFO("Levenshtein distance: " + std::to_string(levenshteinDist) + "/" + std::to_string(targetSequence.length()) + " (ratio: " + std::to_string(levenshteinRatio) + ")");
    }
    
    // Calculate component scores
    double coverageScore = calculateSpectrumCoverageScore(genes, instance);
    double connectivityScore = calculateConnectivityScore(individual, instance);
    double lengthPenalty = calculateLengthPenalty(genes.size(), instance.getN());

    // Use sigmoid normalization for better score differentiation
    auto sigmoid = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };
    
    // Normalize scores to [0,1] range with sigmoid
    double normalizedCoverage = sigmoid((coverageScore - 0.5) * 10);
    double normalizedConnectivity = sigmoid((connectivityScore - 0.5) * 10);
    double normalizedLength = sigmoid((lengthPenalty - 0.5) * 10);

    // Calculate weighted sum with adjusted weights
    double fitness = (normalizedCoverage * 0.5) + 
                    (normalizedConnectivity * 0.3) + 
                    (normalizedLength * 0.2);

    // More lenient validation thresholds
    bool isValid = coverageScore >= 0.7 && 
                  connectivityScore >= 0.6 && 
                  lengthPenalty >= 0.6;

    // Less severe penalty for invalid solutions
    double finalFitness = isValid ? fitness : fitness * 0.8;

    // Add bonus for solutions close to being valid
    if (!isValid && coverageScore >= 0.6 && connectivityScore >= 0.6 && lengthPenalty >= 0.6) {
        finalFitness *= 1.1;  // 10% bonus for almost valid solutions
    }

    LOG_DEBUG("Coverage: " + std::to_string(coverageScore) + 
              ", Connectivity: " + std::to_string(connectivityScore) + 
              ", Length: " + std::to_string(lengthPenalty) + 
              ", Valid: " + std::to_string(isValid));

    return finalFitness;
}

double OptimizedGraphBasedFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    if (!individual) {
        return 0.0;
    }

    const auto& genes = individual->getGenes();
    const auto& spectrum = instance.getSpectrum();
    if (genes.empty()) {
        return 0.0;
    }

    double score = 0.0;
    int k = instance.getK();
    
    // Build adjacency matrix for k-mer overlap graph
    auto adjacencyMatrix = buildAdjacencyMatrix(instance);
    
    // Calculate connectivity based on k-mer overlaps with partial credit
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i + 1] < 0 || genes[i + 1] >= static_cast<int>(spectrum.size())) {
            continue;
        }
        
        const auto& edge = adjacencyMatrix[genes[i]][genes[i + 1]];
        if (edge.exists) {
            score += edge.weight;
        } else {
            // Give partial credit for near-matches
            int partialWeight = calculatePartialOverlapWeight(spectrum[genes[i]], spectrum[genes[i + 1]], k);
            score += partialWeight * 0.5;  // 50% credit for partial overlaps
        }
    }
    
    // Normalize by maximum possible score
    double maxScore = (genes.size() - 1) * k;  // Maximum possible overlap score
    return maxScore > 0.0 ? score / maxScore : 0.0;
}

double OptimizedGraphBasedFitness::calculateSpectrumCoverageScore(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    const auto& spectrum = instance.getSpectrum();
    if (genes.empty()) {
        return 0.0;
    }

    // Track which k-mers from the spectrum have been used
    std::vector<bool> used(spectrum.size(), false);
    int matches = 0;

    // Mark each gene's k-mer as used
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] >= 0 && genes[i] < static_cast<int>(spectrum.size()) && !used[genes[i]]) {
            used[genes[i]] = true;
            matches++;
        }
    }

    // Calculate coverage ratio
    double coverageRatio = static_cast<double>(matches) / spectrum.size();

    // Apply progressive penalties for low coverage
    if (coverageRatio < 0.8) {
        coverageRatio *= 0.8;  // Less severe penalty for < 80% coverage
    } else if (coverageRatio < 0.9) {
        coverageRatio *= 0.9;  // Small penalty for < 90% coverage
    }

    return coverageRatio;
}

double OptimizedGraphBasedFitness::calculateLengthPenalty(
    int actualLength,
    int targetLength) const {
    double diff = std::abs(actualLength - targetLength);
    return std::max(0.0, 1.0 - (diff * diff) / (targetLength * targetLength));
}

std::vector<std::vector<PreprocessedEdge>> OptimizedGraphBasedFitness::buildAdjacencyMatrix(
    const DNAInstance& instance) const {
    const auto& spectrum = instance.getSpectrum();
    std::vector<std::vector<PreprocessedEdge>> matrix(spectrum.size(),
        std::vector<PreprocessedEdge>(spectrum.size()));
    
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], instance.getK());
                matrix[i][j] = PreprocessedEdge(j, weight, weight >= instance.getK() - instance.getDeltaK());
            }
        }
    }
    
    return matrix;
}

int OptimizedGraphBasedFitness::calculateEdgeWeight(
    const std::string& from,
    const std::string& to,
    int k) const {
    if (from.length() < static_cast<size_t>(k - 1) || to.length() < static_cast<size_t>(k)) {
        return 0;
    }
    
    std::string suffix = from.substr(from.length() - (k - 1));
    std::string prefix = to.substr(0, k - 1);
    
    return (suffix == prefix) ? k - 1 : 0;
}

int OptimizedGraphBasedFitness::calculatePartialOverlapWeight(
    const std::string& from,
    const std::string& to,
    int k) const {
    if (from.length() < static_cast<size_t>(k - 1) || to.length() < static_cast<size_t>(k)) {
        return 0;
    }
    
    std::string suffix = from.substr(from.length() - (k - 1));
    std::string prefix = to.substr(0, k - 1);
    
    // Count matching characters
    int matches = 0;
    for (size_t i = 0; i < suffix.length() && i < prefix.length(); ++i) {
        if (suffix[i] == prefix[i]) matches++;
    }
    
    return matches;
}

int OptimizedGraphBasedFitness::calculateLevenshteinDistance(
    const std::string& s1, 
    const std::string& s2) const {
    
    std::vector<std::vector<int>> dp(s1.length() + 1, 
                                    std::vector<int>(s2.length() + 1));
    
    for (size_t i = 0; i <= s1.length(); i++) {
        dp[i][0] = i;
    }
    
    for (size_t j = 0; j <= s2.length(); j++) {
        dp[0][j] = j;
    }
    
    for (size_t i = 1; i <= s1.length(); i++) {
        for (size_t j = 1; j <= s2.length(); j++) {
            if (s1[i-1] == s2[j-1]) {
                dp[i][j] = dp[i-1][j-1];
            } else {
                dp[i][j] = 1 + std::min({dp[i-1][j],      // deletion
                                        dp[i][j-1],        // insertion
                                        dp[i-1][j-1]});    // substitution
            }
        }
    }
    
    return dp[s1.length()][s2.length()];
} 