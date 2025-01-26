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
    const auto& spectrum = instance.getSpectrum();
    
    // Calculate component scores
    double coverageScore = calculateSpectrumCoverageScore(genes, instance);
    double connectivityScore = calculateConnectivityScore(individual, instance);
    double lengthPenalty = calculateLengthPenalty(genes.size(), instance.getN());

    // Use sigmoid normalization for better score differentiation
    auto sigmoid = [](double x) { return 1.0 / (1.0 + std::exp(-x)); };
    
    // Normalize scores to [0,1] range with sigmoid
    // Adjusted thresholds and steepness
    double normalizedCoverage = sigmoid((coverageScore - 0.7) * 12);  // Expect at least 70% coverage
    double normalizedConnectivity = sigmoid((connectivityScore - 0.2) * 15);  // More lenient threshold, steeper curve
    double normalizedLength = sigmoid((lengthPenalty - 0.4) * 8);

    // Calculate weighted sum with increased weight on connectivity
    double fitness = (normalizedCoverage * 0.25) + 
                    (normalizedConnectivity * 0.60) +  // Increased from 0.45 to 0.60
                    (normalizedLength * 0.15);

    // Apply quality penalties
    double qualityMultiplier = 1.0;
    
    // Stronger penalties for poor connectivity
    if (connectivityScore < 0.1) {
        qualityMultiplier *= 0.2;  // Severe penalty for very poor connectivity
    } else if (connectivityScore < 0.2) {
        qualityMultiplier *= 0.5;  // Moderate penalty
    }
    
    // Stronger penalties for low coverage
    if (coverageScore < 0.6) {
        qualityMultiplier *= 0.3;  // Severe penalty for very low coverage
    } else if (coverageScore < 0.7) {
        qualityMultiplier *= 0.6;  // Moderate penalty
    }
    
    if (lengthPenalty < 0.3) {
        qualityMultiplier *= 0.6;
    }
    
    // Add bonus for high-quality solutions
    if (coverageScore >= 0.8 && connectivityScore >= 0.3) {
        qualityMultiplier *= 1.3;  // 30% bonus for good solutions
    }
    
    // Additional bonus for very high connectivity
    if (connectivityScore >= 0.4) {
        qualityMultiplier *= 1.5;  // 50% extra bonus for excellent connectivity
    }

    // Apply quality multiplier
    double finalFitness = fitness * qualityMultiplier;

    // Detailed logging
    LOG_DEBUG("Scores - Coverage: " + std::to_string(coverageScore) + 
              ", Connectivity: " + std::to_string(connectivityScore) + 
              ", Length: " + std::to_string(lengthPenalty) + 
              ", Quality: " + std::to_string(qualityMultiplier));

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
    
    // Track consecutive good overlaps for bonus
    int consecutiveGoodOverlaps = 0;
    double chainBonus = 1.0;
    
    // Calculate connectivity with extreme preference for good overlaps
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i + 1] < 0 || genes[i + 1] >= static_cast<int>(spectrum.size())) {
            consecutiveGoodOverlaps = 0;
            chainBonus = 1.0;
            continue;
        }
        
        // Skip if either k-mer is too short
        if (spectrum[genes[i]].length() < static_cast<size_t>(k - 1) ||
            spectrum[genes[i + 1]].length() < static_cast<size_t>(k)) {
            consecutiveGoodOverlaps = 0;
            chainBonus = 1.0;
            continue;
        }
        
        // Check overlap quality
        std::string suffix = spectrum[genes[i]].substr(spectrum[genes[i]].length() - (k - 1));
        std::string prefix = spectrum[genes[i + 1]].substr(0, k - 1);
        
        int mismatches = 0;
        for (size_t j = 0; j < static_cast<size_t>(k - 1); ++j) {
            if (suffix[j] != prefix[j]) mismatches++;
        }
        
        double overlapScore = 0.0;
        
        // Allow slightly more mismatches than deltaK for partial credit
        int allowedMismatches = instance.getDeltaK() + 1;
        
        if (mismatches <= instance.getDeltaK()) {
            // Perfect overlap within deltaK
            overlapScore = 4.0;  // Increased base score
            consecutiveGoodOverlaps++;
            
            // Exponential bonus for consecutive good overlaps
            chainBonus = 1.0 + (0.2 * consecutiveGoodOverlaps);  // 20% increase per consecutive overlap
            overlapScore *= chainBonus;
            
        } else if (mismatches <= allowedMismatches) {
            // Near-perfect overlap gets partial credit
            overlapScore = 1.0 * (1.0 - static_cast<double>(mismatches) / allowedMismatches);
            consecutiveGoodOverlaps = 0;
            chainBonus = 1.0;
        } else {
            // Poor overlap
            overlapScore = 0.1;  // Small baseline score
            consecutiveGoodOverlaps = 0;
            chainBonus = 1.0;
        }
        
        score += overlapScore;
    }
    
    // Normalize by maximum possible score
    double maxScore = (genes.size() - 1) * 4.0;  // Base max score
    return maxScore > 0.0 ? score / maxScore : 0.0;
}

double OptimizedGraphBasedFitness::calculateSpectrumCoverageScore(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    const auto& spectrum = instance.getSpectrum();
    if (genes.empty()) {
        return 0.0;
    }

    // Track which k-mers have been used
    std::vector<bool> used(spectrum.size(), false);
    int validKmers = 0;
    int usedKmers = 0;
    
    // Count valid k-mers and mark used ones
    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (spectrum[i].length() >= static_cast<size_t>(instance.getK())) {
            validKmers++;
        }
    }
    
    for (int gene : genes) {
        if (gene >= 0 && gene < static_cast<int>(spectrum.size()) && !used[gene]) {
            if (spectrum[gene].length() >= static_cast<size_t>(instance.getK())) {
                used[gene] = true;
                usedKmers++;
            }
        }
    }

    // Calculate coverage ratio
    double coverageRatio = static_cast<double>(usedKmers) / validKmers;
    
    // Progressive penalties for coverage
    if (coverageRatio < 0.6) {
        return coverageRatio * 0.5;  // Severe penalty
    } else if (coverageRatio < 0.7) {
        return coverageRatio * 0.7;  // Moderate penalty
    } else if (coverageRatio < 0.8) {
        return coverageRatio * 0.9;  // Light penalty
    } else {
        return coverageRatio;  // No penalty for good coverage
    }
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