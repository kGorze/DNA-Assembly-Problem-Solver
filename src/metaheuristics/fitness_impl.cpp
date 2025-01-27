#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/dna_utils.h"
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
    
    // Calculate detailed component scores
    double baseConnectivity = calculateConnectivityScore(solution, instance);
    double baseCoverage = calculateSpectrumCoverageScore(dna, instance);
    double lengthPenalty = calculateLengthPenalty(dna.size(), instance.getN());
    
    // Calculate sub-scores for connectivity
    auto [overlapQuality, connectionCount] = calculateDetailedConnectivity(solution, instance);
    
    // Calculate sub-scores for coverage
    auto [exactMatches, partialMatches] = calculateDetailedCoverage(dna, instance);
    
    // Calculate Levenshtein distance score if original DNA is available
    double levenshteinScore = 0.0;
    const std::string& originalDNA = instance.getOriginalDNA();
    if (!originalDNA.empty()) {
        int distance = calculateLevenshteinDistance(dnaStr, originalDNA);
        levenshteinScore = 1.0 - std::min(1.0, static_cast<double>(distance) / originalDNA.length());
    }
    
    // Combine scores with emphasis on incremental improvements
    double fitness = 0.0;
    
    // Connectivity component (40%)
    double connectivityScore = (0.6 * baseConnectivity) +  // Base connectivity
                             (0.3 * overlapQuality) +      // Quality of overlaps
                             (0.1 * connectionCount);      // Number of valid connections
    fitness += 0.4 * connectivityScore;
    
    // Coverage component (30%)
    double coverageScore = (0.7 * baseCoverage) +         // Base coverage
                          (0.2 * exactMatches) +          // Exact k-mer matches
                          (0.1 * partialMatches);         // Partial matches
    fitness += 0.3 * coverageScore;
    
    // Length component (20%)
    fitness += 0.2 * lengthPenalty;
    
    // Levenshtein component (10% if available)
    if (!originalDNA.empty()) {
        fitness += 0.1 * levenshteinScore;
    }
    
    return fitness;
}

double SimpleFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    
    const auto& genes = solution->getGenes();
    if (genes.size() < 2) return 0.0;
    
    int validConnections = 0;
    int totalConnections = genes.size() - 1;
    double totalOverlapQuality = 0.0;
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        const auto& spectrum = instance.getSpectrum();
        if (genes[i] >= 0 && static_cast<size_t>(genes[i]) < spectrum.size() &&
            genes[i+1] >= 0 && static_cast<size_t>(genes[i+1]) < spectrum.size()) {
            
            std::string current = spectrum[genes[i]];
            std::string next = spectrum[genes[i + 1]];
            int k = instance.getK();
            
            if (static_cast<int>(current.length()) >= k-1 && 
                static_cast<int>(next.length()) >= k-1) {
                std::string suffix = current.substr(current.length() - (k-1));
                std::string prefix = next.substr(0, k-1);
                
                // Calculate overlap quality
                int matches = 0;
                for (int j = 0; j < k-1; ++j) {
                    if (suffix[j] == prefix[j]) matches++;
                }
                
                double overlapQuality = static_cast<double>(matches) / (k-1);
                if (overlapQuality >= 0.8) {  // High quality overlap threshold
                    validConnections++;
                    totalOverlapQuality += overlapQuality;
                }
            }
        }
    }
    
    // Combine both connectivity ratio and overlap quality
    double connectivityRatio = totalConnections > 0 ? 
        static_cast<double>(validConnections) / totalConnections : 0.0;
    double avgOverlapQuality = validConnections > 0 ? 
        totalOverlapQuality / validConnections : 0.0;
    
    return 0.7 * connectivityRatio + 0.3 * avgOverlapQuality;
}

double SimpleFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const {
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty() || dna.empty()) return 0.0;  // No minimum score
    
    int covered = 0;
    int k = instance.getK();
    
    // Only count exact matches
    for (size_t i = 0; i <= dna.size() - k; ++i) {
        std::string kmer(dna.begin() + i, dna.begin() + i + k);
        if (std::find(spectrum.begin(), spectrum.end(), kmer) != spectrum.end()) {
            covered++;
        }
    }
    
    return static_cast<double>(covered) / spectrum.size();
}

double SimpleFitness::calculateLengthPenalty(int actualLength, int targetLength) const {
    if (targetLength <= 0) return 0.0;
    
    double ratio = static_cast<double>(actualLength) / targetLength;
    
    // Progressive penalty based on deviation severity
    if (ratio < 0.5 || ratio > 1.5) {
        return 0.1;  // Severe penalty for extreme deviation
    } else if (ratio < 0.7 || ratio > 1.3) {
        return 0.4;  // Significant penalty for large deviation
    } else if (ratio < 0.9 || ratio > 1.1) {
        return 0.7;  // Moderate penalty for medium deviation
    } else {
        // Small deviations get scaled penalty
        return 1.0 - std::abs(1.0 - ratio);
    }
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

ConnectivityMetrics SimpleFitness::calculateDetailedConnectivity(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    
    const auto& genes = solution->getGenes();
    if (genes.size() < 2) return ConnectivityMetrics();
    
    double totalOverlapQuality = 0.0;
    int validConnections = 0;
    std::vector<double> overlapScores;
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        const auto& spectrum = instance.getSpectrum();
        if (genes[i] >= 0 && static_cast<size_t>(genes[i]) < spectrum.size() &&
            genes[i+1] >= 0 && static_cast<size_t>(genes[i+1]) < spectrum.size()) {
            
            std::string current = spectrum[genes[i]];
            std::string next = spectrum[genes[i + 1]];
            int k = instance.getK();
            
            if (static_cast<int>(current.length()) >= k-1 && 
                static_cast<int>(next.length()) >= k-1) {
                std::string suffix = current.substr(current.length() - (k-1));
                std::string prefix = next.substr(0, k-1);
                
                // Calculate detailed overlap quality
                int matches = 0;
                for (int j = 0; j < k-1; ++j) {
                    if (suffix[j] == prefix[j]) matches++;
                }
                
                double overlapQuality = static_cast<double>(matches) / (k-1);
                overlapScores.push_back(overlapQuality);
                
                // Count any non-zero overlap as a valid connection
                if (overlapQuality > 0) {
                    validConnections++;
                    totalOverlapQuality += overlapQuality;
                }
            }
        }
    }
    
    // Calculate normalized metrics
    double avgOverlapQuality = validConnections > 0 ? totalOverlapQuality / validConnections : 0.0;
    double connectionRatio = static_cast<double>(validConnections) / (genes.size() - 1);
    
    return ConnectivityMetrics(avgOverlapQuality, connectionRatio);
}

CoverageMetrics SimpleFitness::calculateDetailedCoverage(
    const std::vector<char>& dna,
    const DNAInstance& instance) const {
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty() || dna.empty()) return CoverageMetrics();
    
    int exactMatches = 0;
    int partialMatches = 0;
    int k = instance.getK();
    
    for (const auto& targetKmer : spectrum) {
        bool foundExact = false;
        bool foundPartial = false;
        
        // Check each possible k-mer in the DNA sequence
        for (size_t i = 0; i + static_cast<size_t>(k) <= dna.size(); ++i) {
            std::string kmer(dna.begin() + i, dna.begin() + i + k);
            
            if (kmer == targetKmer) {
                exactMatches++;
                foundExact = true;
                break;
            }
            
            // Count partial matches (at least 70% matching bases)
            if (!foundPartial) {
                int matches = 0;
                for (int j = 0; j < k; ++j) {
                    if (kmer[static_cast<size_t>(j)] == targetKmer[static_cast<size_t>(j)]) matches++;
                }
                if (static_cast<double>(matches) / k >= 0.7) {
                    partialMatches++;
                    foundPartial = true;
                }
            }
        }
    }
    
    return CoverageMetrics(
        static_cast<double>(exactMatches) / spectrum.size(),
        static_cast<double>(partialMatches) / spectrum.size()
    );
}

// Implementation of OptimizedGraphBasedFitness methods
double OptimizedGraphBasedFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    [[maybe_unused]] std::shared_ptr<IRepresentation> representation) const {
    if (!solution) return 0.0;
    
    // Try to get cached fitness first
    if (auto cache = std::dynamic_pointer_cast<IPopulationCache>(m_config.getCache())) {
        double cachedFitness = cache->getOrCalculateFitness(solution, instance);
        
        // Check for stagnation and adapt parameters if needed
        if (m_stagnationCount > 3) {  // If stagnated for more than 3 generations
            // Increase exploration through mutation rate
            if (auto mutationOp = m_config.getMutation()) {
                mutationOp->setMutationRate(std::min(0.4, mutationOp->getMutationRate() * 1.5));
            }
            
            // Increase replacement ratio for more population turnover
            if (auto replacementOp = m_config.getReplacement()) {
                replacementOp->setReplacementRatio(std::min(0.9, replacementOp->getReplacementRatio() * 1.2));
            }
            
            // Consider random restart if stagnation persists
            if (m_stagnationCount > 10) {
                m_needsRestart = true;
                m_stagnationCount = 0;  // Reset counter after triggering restart
            }
        } else if (m_stagnationCount == 0) {  // If not stagnating, gradually return to base values
            if (auto mutationOp = m_config.getMutation()) {
                mutationOp->setMutationRate(std::max(0.1, mutationOp->getMutationRate() * 0.9));
            }
            if (auto replacementOp = m_config.getReplacement()) {
                replacementOp->setReplacementRatio(std::max(0.7, replacementOp->getReplacementRatio() * 0.95));
            }
        }
        
        return cachedFitness;
    }
    
    // Calculate component scores
    double edgeQualityScore = calculateEdgeQuality(solution, instance);
    double lengthScore = calculateLength(solution, instance);
    
    // Calculate theoretical maximum scores
    double maxEdgeQuality = 1.0;  // Perfect overlaps throughout
    double maxLength = 1.0;       // Exact target length
    
    // Normalize scores relative to theoretical maximums
    edgeQualityScore = std::min(1.0, edgeQualityScore / maxEdgeQuality);
    lengthScore = std::min(1.0, lengthScore / maxLength);
    
    // Apply progressive scaling to reward incremental improvements
    edgeQualityScore = std::pow(edgeQualityScore, 0.8);  // Less punishing for small improvements
    lengthScore = std::pow(lengthScore, 0.9);            // Slightly less punishing for length
    
    // Weight the components with emphasis on edge quality
    double weightedScore = (0.7 * edgeQualityScore + 0.3 * lengthScore);
    
    // Add bonuses for significant achievements
    if (edgeQualityScore > 0.6) {  // Bonus for good edge quality
        weightedScore *= 1.1;  // 10% bonus
    }
    if (edgeQualityScore > 0.8) {  // Extra bonus for excellent edge quality
        weightedScore *= 1.1;  // Additional 10% bonus
    }
    if (lengthScore > 0.8) {  // Bonus for good length
        weightedScore *= 1.05;  // 5% bonus
    }
    
    // Calculate diversity contribution if we have a population
    if (!m_currentPopulation.empty()) {
        double diversityBonus = 0.0;
        double minDistance = 1.0;
        
        // Find minimum distance to any other solution
        for (const auto& other : m_currentPopulation) {
            if (other && other.get() != solution.get()) {
                double distance = calculateDistance(solution, other);
                minDistance = std::min(minDistance, distance);
            }
        }
        
        // Add diversity bonus (up to 10% extra)
        diversityBonus = minDistance * 0.1;
        weightedScore *= (1.0 + diversityBonus);
    }
    
    // Ensure score is properly bounded
    return std::clamp(weightedScore, 0.0, 1.0);
}

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
            
            // Progressive scoring for edge quality
            if (edge.weight > 0) {
                double normalizedWeight = static_cast<double>(edge.weight) / instance.getK();
                // Apply progressive scaling to reward better overlaps more
                double scaledWeight = std::pow(normalizedWeight, 0.8);
                totalQuality += scaledWeight;
                validEdges++;
            } else if (std::abs(genes[i] - genes[i + 1]) <= instance.getDeltaK()) {
                // Small partial credit for close indices
                totalQuality += 0.2;
                validEdges++;
            }
        }
    }
    
    // Calculate final score with minimum threshold
    double score = validEdges > 0 ? totalQuality / (genes.size() - 1) : 0.0;
    return std::max(0.1, score);  // Ensure minimum non-zero score
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
    const int k = instance.getK();
    const int deltaK = instance.getDeltaK();
    
    std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix(
        spectrum.size(),
        std::vector<PreprocessedEdge>(spectrum.size())
    );
    
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
                // More lenient edge validity check
                bool isValid = weight > 0 || std::abs(static_cast<int>(i) - static_cast<int>(j)) <= deltaK + 3;
                adjacencyMatrix[i][j] = PreprocessedEdge(j, weight, isValid);
            }
        }
    }
    
    return adjacencyMatrix;
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

int OptimizedGraphBasedFitness::calculateEdgeWeight(
    const std::string& from, const std::string& to, int k) const {
    if (from.empty() || to.empty()) return 0;
    
    // Try different overlap lengths
    for (int overlap = k - 1; overlap > 0; overlap--) {
        if (static_cast<int>(from.length()) < overlap || 
            static_cast<int>(to.length()) < overlap) continue;
        
        std::string suffix = from.substr(from.length() - overlap);
        std::string prefix = to.substr(0, overlap);
        
        // Count matches
        int matches = 0;
        for (int i = 0; i < overlap; i++) {
            if (suffix[i] == prefix[i]) matches++;
        }
        
        // More lenient scoring
        double matchRatio = static_cast<double>(matches) / overlap;
        if (matchRatio >= 0.8) return k;      // Excellent match
        else if (matchRatio >= 0.6) return k-1;  // Good match
        else if (matchRatio >= 0.4) return k-2;  // Acceptable match
    }
    
    return 1;  // Give minimum score instead of 0
} 