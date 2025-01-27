#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/individual.h"
#include <algorithm>
#include <random>
#include <unordered_set>
#include <limits>
#include <unordered_map>

// Forward declaration of helper function
std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size);

std::vector<std::shared_ptr<Individual>> OrderCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    // Basic validation
    if (parents.size() < 2) {
        LOG_WARNING("OrderCrossover: Not enough parents for crossover");
        return parents;
    }
    
    if (!parents[0] || !parents[1]) {
        LOG_WARNING("OrderCrossover: Null parent received");
        return parents;
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    
    if (parent1Genes.empty() || parent2Genes.empty()) {
        LOG_WARNING("OrderCrossover: Empty genes in parent");
        return parents;
    }
    
    // Create offspring
    std::vector<std::shared_ptr<Individual>> offspring;
    
    // Get size from parent genes
    size_t size = parent1Genes.size();
    
    // Try to create valid offspring multiple times
    const int maxAttempts = 5;
    int attempts = 0;
    bool success = false;
    
    while (attempts < maxAttempts && !success) {
        try {
            // Perform crossover operation with size parameter
            std::vector<int> child1Genes = performOrderCrossover(parent1Genes, parent2Genes, size);
            std::vector<int> child2Genes = performOrderCrossover(parent2Genes, parent1Genes, size);
            
            // Basic validation of gene indices
            bool genesValid = true;
            for (int gene : child1Genes) {
                if (gene < 0 || gene >= static_cast<int>(instance.getSpectrum().size())) {
                    genesValid = false;
                    break;
                }
            }
            
            if (!genesValid) {
                LOG_WARNING("OrderCrossover: Invalid gene indices in offspring");
                attempts++;
                continue;
            }
            
            // Create offspring individuals
            auto child1 = std::make_shared<Individual>(child1Genes);
            auto child2 = std::make_shared<Individual>(child2Genes);
            
            // Validate offspring
            if (representation->isValid(child1, instance) && representation->isValid(child2, instance)) {
                offspring.push_back(child1);
                offspring.push_back(child2);
                success = true;
            } else {
                LOG_WARNING("OrderCrossover: Invalid offspring produced");
                attempts++;
            }
            
        } catch (const std::exception& e) {
            LOG_WARNING("OrderCrossover: Exception during crossover: " + std::string(e.what()));
            attempts++;
        }
    }
    
    // If no valid offspring were produced, return copies of parents
    if (offspring.empty()) {
        LOG_WARNING("OrderCrossover: Failed to produce valid offspring after " + std::to_string(maxAttempts) + " attempts");
        return parents;
    }
    
    return offspring;
}

// EdgeTable implementation
EdgeRecombination::EdgeTable::EdgeTable(const std::vector<std::shared_ptr<Individual>>& parents, const DNAInstance& instance) {
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
                int overlapQuality = calculateEdgeWeight(spectrum[prev], spectrum[current], k);
                if (overlapQuality > 0) {
                    edges[current].push_back({prev, overlapQuality});
                }
            }
            
            if (current >= 0 && static_cast<size_t>(current) < spectrum.size() &&
                next >= 0 && static_cast<size_t>(next) < spectrum.size()) {
                // Calculate overlap quality for current -> next
                int overlapQuality = calculateEdgeWeight(spectrum[current], spectrum[next], k);
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

std::vector<EdgeRecombination::EdgeInfo> EdgeRecombination::EdgeTable::getNeighbors(int node) const {
    auto it = edges.find(node);
    return it != edges.end() ? it->second : std::vector<EdgeInfo>();
}

void EdgeRecombination::EdgeTable::removeNode(int node) {
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

bool EdgeRecombination::EdgeTable::hasNode(int node) const {
    return edges.find(node) != edges.end();
}

std::vector<std::shared_ptr<Individual>> EdgeRecombination::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    [[maybe_unused]] std::shared_ptr<IRepresentation> representation) {
    if (parents.size() < 2) {
        LOG_WARNING("Insufficient parents for crossover");
        return parents;
    }

    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);

    const auto& parent1 = parents[0]->getGenes();
    const auto& parent2 = parents[1]->getGenes();
    
    if (parent1.empty() || parent2.empty()) {
        LOG_WARNING("Empty parent genes in crossover");
        return parents;
    }

    // Build edge table with DNA overlap information
    EdgeTable edgeTable(parents, instance);
    
    // Create two offspring
    for (int k = 0; k < 2; k++) {
        std::vector<int> childGenes;
        childGenes.reserve(parent1.size());
        std::unordered_set<int> used;
        
        // Start with a random gene from first parent
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(0, parent1.size() - 1);
        int current = parent1[dist(gen)];
        childGenes.push_back(current);
        used.insert(current);
        
        // Build the rest of the path prioritizing high-quality overlaps
        while (childGenes.size() < parent1.size()) {
            const auto neighbors = edgeTable.getNeighbors(current);
            bool foundNext = false;
            
            // First try: Find unused neighbor with best overlap quality
            for (const auto& neighbor : neighbors) {
                if (used.find(neighbor.node) == used.end()) {
                    childGenes.push_back(neighbor.node);
                    used.insert(neighbor.node);
                    current = neighbor.node;
                    edgeTable.removeNode(current);
                    foundNext = true;
                    break;
                }
            }
            
            if (!foundNext) {
                // Second try: Use any unused gene from parents, prioritizing those with some overlap
                std::vector<std::pair<int, int>> candidates; // {gene, overlap_quality}
                for (int gene : parent1) {
                    if (used.find(gene) == used.end()) {
                        int quality = calculatePartialOverlapWeight(
                            instance.getSpectrum()[current],
                            instance.getSpectrum()[gene],
                            instance.getK());
                        candidates.push_back({gene, quality});
                    }
                }
                for (int gene : parent2) {
                    if (used.find(gene) == used.end()) {
                        int quality = calculatePartialOverlapWeight(
                            instance.getSpectrum()[current],
                            instance.getSpectrum()[gene],
                            instance.getK());
                        candidates.push_back({gene, quality});
                    }
                }
                
                if (!candidates.empty()) {
                    // Sort by overlap quality
                    std::sort(candidates.begin(), candidates.end(),
                             [](const auto& a, const auto& b) {
                                 return a.second > b.second;
                             });
                    int next = candidates[0].first;
                    childGenes.push_back(next);
                    used.insert(next);
                    current = next;
                    edgeTable.removeNode(current);
                } else {
                    // Last resort: Use any unused value
                    for (size_t i = 0; i < parent1.size(); i++) {
                        if (used.find(i) == used.end()) {
                            childGenes.push_back(i);
                            used.insert(i);
                            current = i;
                            edgeTable.removeNode(current);
                            break;
                        }
                    }
                }
            }
        }
        
        auto child = std::make_shared<Individual>(childGenes);
        offspring.push_back(child);
    }
    
    return offspring;
}

std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dist(0, size - 1);
    
    // Select two random crossover points
    size_t point1 = dist(gen);
    size_t point2 = dist(gen);
    if (point1 > point2) std::swap(point1, point2);
    
    // Create child with same size as parents
    std::vector<int> child(size, -1);
    std::unordered_set<int> used;
    
    // Copy segment between crossover points from parent1
    for (size_t i = point1; i <= point2; ++i) {
        child[i] = parent1[i];
        used.insert(parent1[i]);
    }
    
    // Fill remaining positions with genes from parent2 in order
    size_t curr = (point2 + 1) % size;
    size_t p2pos = (point2 + 1) % size;
    
    // Continue until we reach point1 again
    while (curr != point1) {
        bool foundUnused = false;
        
        // Try one complete cycle through parent2 starting from p2pos
        for (size_t i = 0; i < size; ++i) {
            size_t checkPos = (p2pos + i) % size;
            if (used.find(parent2[checkPos]) == used.end()) {
                child[curr] = parent2[checkPos];
                used.insert(parent2[checkPos]);
                p2pos = (checkPos + 1) % size;
                foundUnused = true;
                break;
            }
        }
        
        // If no unused gene found in parent2, find any unused value
        if (!foundUnused) {
            for (size_t val = 0; val < size; ++val) {
                if (used.find(val) == used.end()) {
                    child[curr] = val;
                    used.insert(val);
                    break;
                }
            }
        }
        
        curr = (curr + 1) % size;
    }
    
    return child;
}

std::vector<std::shared_ptr<Individual>> PMXCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (parents.size() < 2) {
        LOG_ERROR("Not enough parents for crossover");
        return {};
    }
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(parents.size());
    
    for (size_t i = 0; i < parents.size() - 1; i += 2) {
        const auto& parent1 = parents[i];
        const auto& parent2 = parents[i + 1];
        
        if (!parent1 || !parent2) {
            LOG_ERROR("Null parent in crossover");
            continue;
        }
        
        const auto& genes1 = parent1->getGenes();
        const auto& genes2 = parent2->getGenes();
        
        if (genes1.empty() || genes2.empty()) {
            LOG_ERROR("Empty parent in crossover");
            continue;
        }
        
        size_t size = std::min(genes1.size(), genes2.size());
        
        // Generate two random crossover points
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<size_t> dist(0, size - 1);
        size_t point1 = dist(gen);
        size_t point2 = dist(gen);
        if (point1 > point2) std::swap(point1, point2);
        
        // Create offspring
        std::vector<int> child1(size, -1);
        std::vector<int> child2(size, -1);
        
        // Track used values for each child
        std::unordered_set<int> used1, used2;
        
        // Copy the mapping section
        for (size_t j = point1; j <= point2; ++j) {
            child1[j] = genes2[j];
            child2[j] = genes1[j];
            used1.insert(child1[j]);
            used2.insert(child2[j]);
        }
        
        // Create mapping between values in the mapping section
        std::unordered_map<int, int> mapping1, mapping2;
        for (size_t j = point1; j <= point2; ++j) {
            mapping1[genes1[j]] = genes2[j];
            mapping2[genes2[j]] = genes1[j];
        }
        
        // Fill in remaining positions
        for (size_t j = 0; j < size; ++j) {
            if (j >= point1 && j <= point2) continue;
            
            // For child1
            int value1 = genes1[j];
            int safetyCounter1 = 0;
            std::unordered_set<int> visited1;
            
            while (mapping1.find(value1) != mapping1.end() && safetyCounter1 < size) {
                if (visited1.find(value1) != visited1.end()) {
                    // We've found a cycle, break it by using any unused value
                    for (size_t val = 0; val < size; ++val) {
                        if (used1.find(val) == used1.end()) {
                            value1 = val;
                            break;
                        }
                    }
                    break;
                }
                visited1.insert(value1);
                value1 = mapping1[value1];
                safetyCounter1++;
            }
            
            // If we still haven't found a valid value, use any unused value
            if (used1.find(value1) != used1.end()) {
                for (size_t val = 0; val < size; ++val) {
                    if (used1.find(val) == used1.end()) {
                        value1 = val;
                        break;
                    }
                }
            }
            
            child1[j] = value1;
            used1.insert(value1);
            
            // For child2
            int value2 = genes2[j];
            int safetyCounter2 = 0;
            std::unordered_set<int> visited2;
            
            while (mapping2.find(value2) != mapping2.end() && safetyCounter2 < size) {
                if (visited2.find(value2) != visited2.end()) {
                    // We've found a cycle, break it by using any unused value
                    for (size_t val = 0; val < size; ++val) {
                        if (used2.find(val) == used2.end()) {
                            value2 = val;
                            break;
                        }
                    }
                    break;
                }
                visited2.insert(value2);
                value2 = mapping2[value2];
                safetyCounter2++;
            }
            
            // If we still haven't found a valid value, use any unused value
            if (used2.find(value2) != used2.end()) {
                for (size_t val = 0; val < size; ++val) {
                    if (used2.find(val) == used2.end()) {
                        value2 = val;
                        break;
                    }
                }
            }
            
            child2[j] = value2;
            used2.insert(value2);
        }
        
        auto offspring1 = std::make_shared<Individual>(child1);
        auto offspring2 = std::make_shared<Individual>(child2);
        
        if (representation->isValid(offspring1, instance)) {
            offspring.push_back(offspring1);
        }
        if (representation->isValid(offspring2, instance)) {
            offspring.push_back(offspring2);
        }
    }
    
    // Handle odd number of parents
    if (parents.size() % 2 == 1) {
        offspring.push_back(parents.back());
    }
    
    return offspring;
}

std::vector<DNAAlignmentCrossover::AlignmentSegment> 
DNAAlignmentCrossover::findAlignmentSegments(
    const std::vector<int>& genes,
    const DNAInstance& instance,
    int minLength) const {
    
    std::vector<AlignmentSegment> segments;
    const auto& spectrum = instance.getSpectrum();
    
    // Find all possible segments of length >= minLength
    for (size_t start = 0; start < genes.size(); start++) {
        for (size_t len = minLength; start + len <= genes.size(); len++) {
            double quality = calculateSegmentQuality(genes, start, len, instance);
            if (quality > 0.0) {  // Only keep segments with good overlaps
                std::vector<int> segmentGenes(genes.begin() + start, 
                                            genes.begin() + start + len);
                segments.emplace_back(std::move(segmentGenes), quality);
            }
        }
    }
    
    // Sort segments by quality
    std::sort(segments.begin(), segments.end());
    return segments;
}

double DNAAlignmentCrossover::calculateSegmentQuality(
    const std::vector<int>& genes,
    size_t start,
    size_t length,
    const DNAInstance& instance) const {
    
    if (start + length > genes.size()) return 0.0;
    
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    double totalQuality = 0.0;
    int validOverlaps = 0;
    
    // Check overlaps between consecutive k-mers
    for (size_t i = start; i < start + length - 1; i++) {
        if (genes[i] >= 0 && static_cast<size_t>(genes[i]) < spectrum.size() &&
            genes[i+1] >= 0 && static_cast<size_t>(genes[i+1]) < spectrum.size()) {
            
            int overlapQuality = calculateEdgeWeight(
                spectrum[genes[i]], spectrum[genes[i+1]], k);
            
            if (overlapQuality > 0) {
                totalQuality += static_cast<double>(overlapQuality) / k;
                validOverlaps++;
            }
        }
    }
    
    // Return normalized quality score
    if (validOverlaps == 0) return 0.0;
    return (totalQuality / validOverlaps) * (static_cast<double>(validOverlaps) / (length - 1));
}

std::vector<int> DNAAlignmentCrossover::mergeSegments(
    const std::vector<AlignmentSegment>& segments1,
    const std::vector<AlignmentSegment>& segments2,
    size_t targetLength) const {
    
    std::vector<int> merged;
    merged.reserve(targetLength);
    std::unordered_set<int> used;
    
    // Try to use high-quality segments alternately from both parents
    size_t i1 = 0, i2 = 0;
    bool useFirst = true;
    
    while (merged.size() < targetLength && (i1 < segments1.size() || i2 < segments2.size())) {
        const auto& currentSegments = useFirst ? segments1 : segments2;
        size_t& currentIndex = useFirst ? i1 : i2;
        
        // Find next valid segment
        while (currentIndex < currentSegments.size()) {
            const auto& segment = currentSegments[currentIndex];
            bool canUseSegment = true;
            
            // Check if segment genes are already used
            for (int gene : segment.genes) {
                if (used.count(gene) > 0) {
                    canUseSegment = false;
                    break;
                }
            }
            
            if (canUseSegment) {
                // Add segment genes
                for (int gene : segment.genes) {
                    if (merged.size() < targetLength) {
                        merged.push_back(gene);
                        used.insert(gene);
                    }
                }
                break;
            }
            currentIndex++;
        }
        
        useFirst = !useFirst;  // Alternate between parents
    }
    
    // Fill remaining positions with unused genes from parents
    if (merged.size() < targetLength) {
        std::vector<int> unusedGenes;
        for (const auto& segments : {segments1, segments2}) {
            for (const auto& segment : segments) {
                for (int gene : segment.genes) {
                    if (used.count(gene) == 0) {
                        unusedGenes.push_back(gene);
                        used.insert(gene);
                    }
                }
            }
        }
        
        // Add unused genes until we reach target length
        size_t i = 0;
        while (merged.size() < targetLength && i < unusedGenes.size()) {
            merged.push_back(unusedGenes[i++]);
        }
    }
    
    return merged;
}

std::vector<std::shared_ptr<Individual>> DNAAlignmentCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !representation) {
        LOG_WARNING("Invalid input for DNA alignment crossover");
        return parents;
    }
    
    try {
        auto [genes1, genes2] = validateAndGetGenes(parents);
        
        // Find alignment segments in both parents
        auto segments1 = findAlignmentSegments(genes1, instance);
        auto segments2 = findAlignmentSegments(genes2, instance);
        
        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);
        
        // Create two offspring with different segment combinations
        for (int i = 0; i < 2; i++) {
            auto childGenes = mergeSegments(segments1, segments2, genes1.size());
            auto child = std::make_shared<Individual>(std::move(childGenes));
            offspring.push_back(child);
            
            // Swap segment orders for second child
            std::swap(segments1, segments2);
        }
        
        return offspring;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in DNA alignment crossover: " + std::string(e.what()));
        return parents;
    }
}

// Helper functions for DNA overlap calculations
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