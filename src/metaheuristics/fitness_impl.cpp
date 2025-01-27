#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <unordered_map>

// Helper functions for edge calculations
int calculateEdgeWeight(const std::string& from, const std::string& to, int k) {
    if (from.empty() || to.empty()) return 0;
    
    // Handle variable length k-mers
    int minLength = std::min(from.length(), to.length());
    int maxOverlap = std::min(minLength - 1, k - 1);
    
    if (maxOverlap <= 0) return 0;
    
    std::string suffix = from.substr(from.length() - maxOverlap);
    std::string prefix = to.substr(0, maxOverlap);
    
    // Count matches and calculate quality score
    int matches = 0;
    for (int i = 0; i < maxOverlap; i++) {
        if (suffix[i] == prefix[i]) matches++;
    }
    
    // Return weighted score based on match quality
    if (matches == maxOverlap) {
        return k;  // Perfect match gets full score
    } else if (matches >= maxOverlap - 1) {
        return k - 1;  // One mismatch gets high score
    } else if (matches >= maxOverlap - 2) {
        return k - 2;  // Two mismatches gets medium score
    } else if (matches >= maxOverlap / 2) {
        return k - 3;  // Partial match gets low score
    }
    return 0;  // Poor match gets no score
}

int calculatePartialOverlapWeight(const std::string& from, const std::string& to, int k) {
    if (from.empty() || to.empty()) return 0;
    
    // Handle variable length k-mers
    int minLength = std::min(from.length(), to.length());
    int maxOverlap = std::min(minLength - 1, k - 1);
    
    if (maxOverlap <= 0) return 0;
    
    std::string suffix = from.substr(from.length() - maxOverlap);
    std::string prefix = to.substr(0, maxOverlap);
    
    // Count matches and calculate quality score
    int matches = 0;
    int mismatches = 0;
    for (int i = 0; i < maxOverlap; i++) {
        if (suffix[i] == prefix[i]) {
            matches++;
        } else {
            mismatches++;
        }
    }
    
    // Calculate overlap quality score
    double matchRatio = static_cast<double>(matches) / maxOverlap;
    if (matchRatio >= 0.9) {
        return maxOverlap;  // Excellent overlap
    } else if (matchRatio >= 0.8) {
        return static_cast<int>(maxOverlap * 0.8);  // Good overlap
    } else if (matchRatio >= 0.7) {
        return static_cast<int>(maxOverlap * 0.6);  // Acceptable overlap
    } else if (matchRatio >= 0.5) {
        return static_cast<int>(maxOverlap * 0.4);  // Poor but potentially useful overlap
    }
    return 0;  // Too many mismatches
}

// Implementation of SimpleFitness methods
double SimpleFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const {
    
    if (!solution || !representation) return 0.0;
    
    // Convert genes to DNA sequence for coverage calculation
    std::string dnaStr = representation->toDNA(solution, instance);
    std::vector<char> dna(dnaStr.begin(), dnaStr.end());
    
    // Calculate component scores with more lenient criteria
    double coverage = calculateSpectrumCoverageScore(dna, instance);
    double connectivity = calculateConnectivityScore(solution, instance);
    
    // Weight coverage and connectivity equally
    double weightedScore = (coverage * 0.5) + (connectivity * 0.5);
    
    // Add bonus for solutions that achieve minimum thresholds
    if (coverage > 0.3 && connectivity > 0.2) {  // Lower thresholds
        weightedScore *= 1.2;  // 20% bonus for balanced solutions
    }
    
    return weightedScore;
}

double SimpleFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    
    if (!solution) return 0.0;
    
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;
    
    int validConnections = 0;
    int totalConnections = genes.size() - 1;
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        // More lenient connectivity check
        if (std::abs(genes[i] - genes[i + 1]) <= instance.getDeltaK() + 2) {  // Allow more distance
            validConnections++;
        }
    }
    
    return totalConnections > 0 ? 
           static_cast<double>(validConnections) / totalConnections : 0.0;
}

double SimpleFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const {
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty() || dna.empty()) return 0.0;
    
    int covered = 0;
    int partialMatches = 0;
    
    for (const auto& kmer : spectrum) {
        // Check for exact matches
        if (std::search(dna.begin(), dna.end(), kmer.begin(), kmer.end()) != dna.end()) {
            covered++;
            continue;
        }
        
        // Check for partial matches (at least 60% similarity)
        for (auto it = dna.begin(); it + kmer.size() <= dna.end(); ++it) {
            int matches = 0;
            for (size_t i = 0; i < kmer.size(); ++i) {
                if (*(it + i) == kmer[i]) matches++;
            }
            if (static_cast<double>(matches) / kmer.size() >= 0.6) {  // More lenient partial match threshold
                partialMatches++;
                break;
            }
        }
    }
    
    // Calculate coverage ratio with partial credit and minimum threshold
    double coverageRatio = (covered + (partialMatches * 0.6)) / spectrum.size();  // More credit for partial matches
    return std::max(0.1, coverageRatio);  // Ensure minimum non-zero fitness
}

double SimpleFitness::calculateLengthPenalty(int actualLength, int targetLength) const {
    double diff = std::abs(actualLength - targetLength);
    return std::max(0.0, 1.0 - diff / targetLength);
}

int SimpleFitness::calculateLevenshteinDistance(
    const std::string& s1, const std::string& s2) const {
    
    std::vector<std::vector<int>> dp(s1.length() + 1, 
                                    std::vector<int>(s2.length() + 1));
    
    for (size_t i = 0; i <= s1.length(); i++) dp[i][0] = i;
    for (size_t j = 0; j <= s2.length(); j++) dp[0][j] = j;
    
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

// Implementation of OptimizedGraphBasedFitness methods
double OptimizedGraphBasedFitness::calculateEdgeQuality(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    if (!individual) return 0.0;
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return 0.0;
    
    auto adjacencyMatrix = buildAdjacencyMatrix(instance);
    double totalQuality = 0.0;
    int validEdges = 0;
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (static_cast<size_t>(genes[i]) < adjacencyMatrix.size() && 
            static_cast<size_t>(genes[i + 1]) < adjacencyMatrix[genes[i]].size()) {
            const auto& edge = adjacencyMatrix[genes[i]][genes[i + 1]];
            if (edge.valid) {
                totalQuality += static_cast<double>(edge.weight) / instance.getK();
                validEdges++;
            }
        }
    }
    
    return validEdges > 0 ? totalQuality / validEdges : 0.0;
}

double OptimizedGraphBasedFitness::calculateLength(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    if (!individual) return 0.0;
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return 0.0;
    
    int totalLength = 0;
    for (size_t i = 0; i < genes.size(); ++i) {
        if (static_cast<size_t>(genes[i]) < instance.getSpectrum().size()) {
            totalLength += instance.getSpectrum()[genes[i]].length();
        }
    }
    
    int totalOverlap = 0;
    auto adjacencyMatrix = buildAdjacencyMatrix(instance);
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (static_cast<size_t>(genes[i]) < adjacencyMatrix.size() && 
            static_cast<size_t>(genes[i + 1]) < adjacencyMatrix[genes[i]].size()) {
            const auto& edge = adjacencyMatrix[genes[i]][genes[i + 1]];
            totalOverlap += edge.overlap;
        }
    }
    
    int finalLength = totalLength - totalOverlap;
    int targetLength = instance.getK() * genes.size(); // Estimate target length based on k-mer size
    
    return calculateLengthPenalty(finalLength, targetLength);
}

double OptimizedGraphBasedFitness::calculateCoverage(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    
    if (!solution) return 0.0;
    
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;
    
    // Calculate coverage based on spectrum matches
    std::vector<char> dna;  // Convert genes to DNA sequence
    for (const auto& gene : genes) {
        if (static_cast<size_t>(gene) < instance.getSpectrum().size()) {
            const auto& kmer = instance.getSpectrum()[gene];
            dna.insert(dna.end(), kmer.begin(), kmer.end());
        }
    }
    
    return calculateSpectrumCoverageScore(dna, instance);
}

double OptimizedGraphBasedFitness::calculateConnectivity(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    
    if (!solution) return 0.0;
    
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;
    
    // Build adjacency matrix for efficient edge weight calculation
    auto adjacencyMatrix = buildAdjacencyMatrix(instance);
    
    double totalWeight = 0.0;
    int validEdges = 0;
    
    // Calculate connectivity based on edge weights
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (static_cast<size_t>(genes[i]) < adjacencyMatrix.size() && 
            static_cast<size_t>(genes[i + 1]) < adjacencyMatrix[genes[i]].size()) {
            const auto& edge = adjacencyMatrix[genes[i]][genes[i + 1]];
            if (edge.valid) {
                totalWeight += edge.weight;
                validEdges++;
            }
        }
    }
    
    return validEdges > 0 ? totalWeight / (validEdges * instance.getK()) : 0.0;
}

std::vector<std::vector<PreprocessedEdge>> OptimizedGraphBasedFitness::buildAdjacencyMatrix(
    const DNAInstance& instance) const {
    
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    // Initialize adjacency matrix
    std::vector<std::vector<PreprocessedEdge>> matrix(
        spectrum.size(),
        std::vector<PreprocessedEdge>(spectrum.size())
    );
    
    // Build edges between k-mers
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {  // Don't create self-loops
                const auto& from = spectrum[i];
                const auto& to = spectrum[j];
                
                int weight = calculateEdgeWeight(from, to, k);
                int overlap = calculatePartialOverlapWeight(from, to, k);
                
                matrix[i][j] = PreprocessedEdge(
                    weight,
                    overlap,
                    weight > 0 || overlap > 0
                );
            }
        }
    }
    
    return matrix;
}

double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const {
    if (!solution) return 0.0;
    
    // Try to get cached fitness first
    if (auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache())) {
        return cache->getOrCalculateFitness(solution, instance);
    }
    
    // Calculate coverage and connectivity scores
    double coverageScore = calculateCoverage(solution, instance);
    double connectivityScore = calculateConnectivity(solution, instance);
    
    // Get the DNA sequence
    std::string dna = representation->toDNA(solution, instance);
    
    // Calculate length penalty - much more forgiving
    double lengthPenalty = 1.0;
    int targetLength = instance.getN();
    if (static_cast<int>(dna.length()) != targetLength) {
        // Linear penalty that's gentler for lengths within ±30% of target
        double lengthDiff = std::abs(static_cast<double>(dna.length()) - targetLength);
        double relativeDiff = lengthDiff / targetLength;
        lengthPenalty = std::max(0.3, 1.0 - (relativeDiff * 0.5));  // More lenient penalty
    }
    
    // Calculate N-base penalty - more forgiving
    double nBasePenalty = 1.0;
    int nCount = std::count(dna.begin(), dna.end(), 'N');
    if (nCount > 0) {
        // Allow up to 20% N bases with minimal penalty
        double nRatio = static_cast<double>(nCount) / dna.length();
        nBasePenalty = std::max(0.4, 1.0 - (nRatio * 0.8));  // More lenient penalty
    }
    
    // Weight the components - emphasize coverage more
    double weightedScore = (0.8 * coverageScore + 0.2 * connectivityScore);
    
    // Apply penalties more gently
    weightedScore *= std::sqrt(lengthPenalty);  // Square root makes penalty more gradual
    weightedScore *= std::sqrt(nBasePenalty);   // Square root makes penalty more gradual
    
    // Add bonuses for partial achievements
    if (coverageScore > 0.1) weightedScore *= 1.1;  // 10% bonus for basic coverage
    if (connectivityScore > 0.1) weightedScore *= 1.1;  // 10% bonus for basic connectivity
    if (coverageScore > 0.3) weightedScore *= 1.2;  // 20% additional bonus for good coverage
    if (connectivityScore > 0.3) weightedScore *= 1.15;  // 15% additional bonus for good connectivity
    
    // Bonus for balanced solutions - more lenient thresholds
    if (coverageScore > 0.2 && connectivityScore > 0.2) {
        weightedScore *= 1.3;  // 30% bonus for balanced solutions
    }
    
    return weightedScore;
}

double OptimizedGraphBasedFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    
    if (!individual) return 0.0;
    
    const auto& genes = individual->getGenes();
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    const int deltaK = instance.getDeltaK();
    const int lNeg = instance.getLNeg();
    const int lPoz = instance.getLPoz();
    
    if (genes.empty() || spectrum.empty()) return 0.0;
    
    double totalScore = 0.0;
    int validEdges = 0;
    int consecutiveGoodOverlaps = 0;  // Track consecutive good overlaps for bonus
    
    // Calculate allowed mismatches based on instance parameters
    int allowedMismatches = deltaK;
    if (lNeg > 0 || lPoz > 0) {
        allowedMismatches += 1;  // Allow one extra mismatch with errors present
    }
    
    for (size_t i = 0; i < genes.size() - 1; i++) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i+1] < 0 || genes[i+1] >= static_cast<int>(spectrum.size())) {
            continue;
        }
        
        const std::string& current = spectrum[genes[i]];
        const std::string& next = spectrum[genes[i+1]];
        
        // Handle variable length k-mers
        int overlapLength = k - 1;
        if (deltaK > 0) {
            overlapLength = std::min({
                k - 1,
                static_cast<int>(current.length()) - 1,
                static_cast<int>(next.length()) - 1
            });
        }
        
        if (static_cast<size_t>(overlapLength + 1) > current.length() || 
            static_cast<size_t>(overlapLength + 1) > next.length()) {
            continue;
        }
        
        // Count mismatches in overlap region
        int mismatches = 0;
        std::string suffix = current.substr(current.length() - overlapLength);
        std::string prefix = next.substr(0, overlapLength);
        
        for (int j = 0; j < overlapLength; j++) {
            if (suffix[j] != prefix[j]) mismatches++;
        }
        
        // Score based on instance parameters
        double edgeScore = 0.0;
        if (mismatches == 0) {
            edgeScore = 1.0;
            consecutiveGoodOverlaps++;
        } else if (mismatches <= allowedMismatches) {
            edgeScore = 1.0 - (static_cast<double>(mismatches) / (allowedMismatches + 1));
            if (edgeScore >= 0.7) consecutiveGoodOverlaps++;
            else consecutiveGoodOverlaps = 0;
        } else {
            edgeScore = 0.1;  // Small score for connected but mismatched
            consecutiveGoodOverlaps = 0;
        }
        
        // Apply bonus for consecutive good overlaps
        if (consecutiveGoodOverlaps >= 3) {
            edgeScore *= (1.0 + 0.1 * std::min(consecutiveGoodOverlaps - 2, 5));
        }
        
        totalScore += edgeScore;
        validEdges++;
    }
    
    // Normalize score
    double normalizedScore = validEdges > 0 ? totalScore / validEdges : 0.0;
    
    // Apply scaling based on instance parameters
    if (deltaK > 0) {
        // More lenient scoring with variable lengths
        normalizedScore = std::pow(normalizedScore, 0.9);
    }
    if (lNeg > 0 || lPoz > 0) {
        // More lenient scoring with errors
        normalizedScore = std::pow(normalizedScore, 0.85);
    }
    
    return normalizedScore;
}

double OptimizedGraphBasedFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const {
    
    const auto& spectrum = instance.getSpectrum();
    int covered = 0;
    
    for (const auto& kmer : spectrum) {
        if (std::search(dna.begin(), dna.end(), kmer.begin(), kmer.end()) != dna.end()) {
            covered++;
        }
    }
    
    return spectrum.empty() ? 0.0 : static_cast<double>(covered) / spectrum.size();
}

double OptimizedGraphBasedFitness::calculateLengthPenalty(
    int actualLength,
    int targetLength) const {
    double diff = std::abs(actualLength - targetLength);
    return std::max(0.0, 1.0 - (diff * diff) / (targetLength * targetLength));
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