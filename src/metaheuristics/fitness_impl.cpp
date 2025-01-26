#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <unordered_map>

double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const {
    if (!individual || !representation) {
        return 0.0;
    }

    const auto& genes = individual->getGenes();
    const auto& spectrum = instance.getSpectrum();
    
    // Calculate component scores with weighted edge approach
    double coverageScore = calculateSpectrumCoverageScore(genes, instance);
    double connectivityScore = calculateConnectivityScore(individual, instance);
    double lengthPenalty = calculateLengthPenalty(genes.size(), instance.getN());

    // Track repeated k-mers
    std::unordered_map<int, int> kmerUsage;
    for (int gene : genes) {
        kmerUsage[gene]++;
    }
    
    // Calculate repetition penalty
    double repetitionPenalty = 0.0;
    for (const auto& [_, count] : kmerUsage) {
        if (count > 1) {
            repetitionPenalty += (count - 1) * 0.1; // Penalty for each repeat
        }
    }

    // Calculate edge weight quality
    double edgeQuality = 0.0;
    int totalEdges = 0;
    const int k = instance.getK();
    const int deltaK = instance.getDeltaK();
    
    for (size_t i = 0; i < genes.size() - 1; i++) {
        if (genes[i] >= 0 && genes[i] < static_cast<int>(spectrum.size()) &&
            genes[i+1] >= 0 && genes[i+1] < static_cast<int>(spectrum.size())) {
            
            const std::string& current = spectrum[genes[i]];
            const std::string& next = spectrum[genes[i+1]];
            
            if (current.length() >= k-1 && next.length() >= k-1) {
                int mismatches = 0;
                std::string suffix = current.substr(current.length() - (k-1));
                std::string prefix = next.substr(0, k-1);
                
                for (int j = 0; j < k-1; j++) {
                    if (suffix[j] != prefix[j]) mismatches++;
                }
                
                // Score edges based on mismatch count
                if (mismatches == 0) {
                    edgeQuality += 1.0; // Perfect match
                } else if (mismatches <= deltaK) {
                    edgeQuality += 0.7; // Good match within deltaK
                } else if (mismatches <= deltaK + 1) {
                    edgeQuality += 0.3; // Acceptable match
                } else {
                    edgeQuality += 0.1; // Poor match but still connected
                }
                totalEdges++;
            }
        }
    }
    
    // Normalize edge quality
    edgeQuality = totalEdges > 0 ? edgeQuality / totalEdges : 0.0;

    // Combine scores with adjusted weights
    double fitness = (coverageScore * 0.3) +           // Coverage of spectrum
                    (connectivityScore * 0.3) +        // Basic connectivity
                    (edgeQuality * 0.25) +             // Quality of connections
                    (lengthPenalty * 0.15);           // Length penalty
                    
    // Apply penalties
    fitness *= (1.0 - std::min(0.5, repetitionPenalty)); // Cap penalty at 50%

    // Only log significant improvements
    static double lastLoggedFitness = 0.0;
    if (fitness > lastLoggedFitness + 0.01) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(4);
        ss << "Improved Solution - Fitness: " << fitness;
        ss << " [Coverage:" << coverageScore;
        ss << " Connectivity:" << connectivityScore;
        ss << " EdgeQuality:" << edgeQuality;
        ss << " Length:" << lengthPenalty;
        if (repetitionPenalty > 0) {
            ss << " RepeatPenalty:" << repetitionPenalty;
        }
        ss << "]";
        LOG_DEBUG(ss.str());
        lastLoggedFitness = fitness;
    }

    return fitness;
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
    const int k = instance.getK();
    const int deltaK = instance.getDeltaK();
    
    // Calculate additional allowed mismatches based on error rates
    const int totalErrors = instance.getLNeg() + instance.getLPoz();
    const int extraMismatches = (deltaK == 0 && totalErrors > 0) ? 
        std::min(2, totalErrors / 10) : 0;
    const int allowedMismatches = deltaK + extraMismatches;
    
    // Track consecutive good overlaps for bonus
    int consecutiveGoodOverlaps = 0;
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i + 1] < 0 || genes[i + 1] >= static_cast<int>(spectrum.size())) {
            continue;
        }

        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i + 1]];

        if (current.length() >= k - 1 && next.length() >= k - 1) {
            std::string suffix = current.substr(current.length() - (k - 1));
            std::string prefix = next.substr(0, k - 1);
            
            int mismatches = 0;
            for (int j = 0; j < k - 1; j++) {
                if (suffix[j] != prefix[j]) mismatches++;
            }
            
            double overlapScore = 0.0;
            
            if (mismatches == 0) {
                // Perfect overlap
                overlapScore = 1.0;
                consecutiveGoodOverlaps++;
                // Bonus for consecutive good overlaps
                overlapScore *= (1.0 + (0.1 * consecutiveGoodOverlaps));
            } else if (mismatches <= deltaK) {
                // Good overlap within deltaK
                overlapScore = 0.8 * (1.0 - static_cast<double>(mismatches) / deltaK);
                consecutiveGoodOverlaps = 0;
            } else if (mismatches <= allowedMismatches) {
                // Acceptable overlap with extra mismatches
                overlapScore = 0.4 * (1.0 - static_cast<double>(mismatches - deltaK) / extraMismatches);
                consecutiveGoodOverlaps = 0;
            } else {
                // Poor overlap but still connected
                overlapScore = 0.1;
                consecutiveGoodOverlaps = 0;
            }
            
            score += overlapScore;
        }
    }
    
    // Normalize by maximum possible score (considering potential bonuses)
    double maxPossibleScore = genes.size() - 1;  // Base score for perfect overlaps
    maxPossibleScore *= 1.5;  // Account for potential consecutive bonuses
    
    return maxPossibleScore > 0.0 ? score / maxPossibleScore : 0.0;
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