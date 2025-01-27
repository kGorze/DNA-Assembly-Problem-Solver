#pragma once

#include <string>
#include <algorithm>

namespace dna_utils {

inline int calculateEdgeWeight(const std::string& from, const std::string& to, int k) {
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

inline int calculatePartialOverlapWeight(const std::string& from, const std::string& to, int k) {
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

} // namespace dna_utils 