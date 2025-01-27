#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/metaheuristics/dna_utils.h"
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
        LOG_WARNING("OrderCrossover: Empty parent genes");
        return parents;
    }
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);
    
    int attempts = 0;
    const int MAX_ATTEMPTS = 3;
    bool success = false;
    
    while (!success && attempts < MAX_ATTEMPTS) {
        try {
            // Create child genes
            auto child1Genes = performOrderCrossover(parent1Genes, parent2Genes, parent1Genes.size());
            auto child2Genes = performOrderCrossover(parent2Genes, parent1Genes, parent2Genes.size());
            
            if (child1Genes.empty() || child2Genes.empty()) {
                LOG_WARNING("OrderCrossover: Invalid gene indices in offspring");
                attempts++;
                continue;
            }
            
            // Create offspring individuals
            auto child1 = std::make_shared<Individual>(child1Genes);
            auto child2 = std::make_shared<Individual>(child2Genes);
            
            // Validate offspring
            bool child1Valid = !child1->getGenes().empty() && representation->isValid(child1, instance);
            bool child2Valid = !child2->getGenes().empty() && representation->isValid(child2, instance);
            
            if (child1Valid) {
                offspring.push_back(child1);
            }
            if (child2Valid) {
                offspring.push_back(child2);
            }
            
            if (!offspring.empty()) {
                success = true;
            } else {
                LOG_WARNING("OrderCrossover: No valid offspring produced");
                attempts++;
            }
            
        } catch (const std::exception& e) {
            LOG_WARNING("OrderCrossover: Exception during crossover: " + std::string(e.what()));
            attempts++;
        }
    }
    
    // If no valid offspring were produced, return copies of parents
    if (offspring.empty()) {
        LOG_WARNING("OrderCrossover: Failed to create valid offspring, returning parents");
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
                int overlapQuality = dna_utils::calculateEdgeWeight(spectrum[prev], spectrum[current], k);
                if (overlapQuality > 0) {
                    edges[current].push_back({prev, overlapQuality});
                }
            }
            
            if (current >= 0 && static_cast<size_t>(current) < spectrum.size() &&
                next >= 0 && static_cast<size_t>(next) < spectrum.size()) {
                // Calculate overlap quality for current -> next
                int overlapQuality = dna_utils::calculateEdgeWeight(spectrum[current], spectrum[next], k);
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
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !parents[0] || !parents[1]) {
        LOG_WARNING("EdgeRecombination: Invalid parents");
        return parents;
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    
    if (parent1Genes.empty() || parent2Genes.empty()) {
        LOG_WARNING("EdgeRecombination: Empty parent genes");
        return parents;
    }
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);
    
    int maxAttempts = 3;
    while (offspring.empty() && maxAttempts > 0) {
        try {
            EdgeTable edgeTable(parents, instance);
            auto& rng = Random::instance();
            
            // Try to create two offspring
            for (int k = 0; k < 2 && offspring.size() < 2; k++) {
                std::vector<int> childGenes;
                childGenes.reserve(parent1Genes.size());
                
                // Start with a random gene from first parent
                int current = parent1Genes[rng.getRandomInt(0, static_cast<int>(parent1Genes.size() - 1))];
                childGenes.push_back(current);
                
                // Build the rest of the path prioritizing high-quality overlaps
                bool validPath = true;
                while (childGenes.size() < parent1Genes.size() && validPath) {
                    const auto neighbors = edgeTable.getNeighbors(current);
                    bool foundNext = false;
                    
                    // First try: Find unused neighbor with best overlap quality
                    for (const auto& neighbor : neighbors) {
                        if (std::find(childGenes.begin(), childGenes.end(), neighbor.node) == childGenes.end()) {
                            childGenes.push_back(neighbor.node);
                            current = neighbor.node;
                            foundNext = true;
                            break;
                        }
                    }
                    
                    if (!foundNext) {
                        // Second try: Use any unused gene from parents
                        std::vector<int> unusedGenes;
                        for (int gene : parent1Genes) {
                            if (std::find(childGenes.begin(), childGenes.end(), gene) == childGenes.end()) {
                                unusedGenes.push_back(gene);
                            }
                        }
                        for (int gene : parent2Genes) {
                            if (std::find(childGenes.begin(), childGenes.end(), gene) == childGenes.end()) {
                                unusedGenes.push_back(gene);
                            }
                        }
                        
                        if (!unusedGenes.empty()) {
                            // Choose random unused gene
                            int idx = rng.getRandomInt(0, static_cast<int>(unusedGenes.size() - 1));
                            childGenes.push_back(unusedGenes[idx]);
                            current = unusedGenes[idx];
                        } else {
                            validPath = false;
                        }
                    }
                }
                
                if (validPath) {
                    auto child = std::make_shared<Individual>(childGenes);
                    if (!child->getGenes().empty() && representation->isValid(child, instance)) {
                        offspring.push_back(child);
                    }
                }
            }
        } catch (const std::exception& e) {
            LOG_WARNING("EdgeRecombination: Exception during crossover: " + std::string(e.what()));
        }
        maxAttempts--;
    }
    
    // If we couldn't create any valid offspring, return copies of parents
    if (offspring.empty()) {
        LOG_WARNING("EdgeRecombination: Failed to create valid offspring, returning parents");
        return parents;
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
    
    // Create child with same size as parents, initialized with first valid value
    std::vector<int> child(size, 0);  // Initialize with 0 instead of -1
    std::unordered_set<int> used;
    
    // Copy segment between crossover points from parent1
    for (size_t i = point1; i <= point2; ++i) {
        child[i] = parent1[i];
        used.insert(parent1[i]);
    }
    
    // Fill remaining positions with genes from parent2 in order
    size_t curr = (point2 + 1) % size;
    size_t p2pos = (point2 + 1) % size;
    
    while (curr != point1) {
        // Try to find an unused gene from parent2
        bool foundUnused = false;
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
        
        // If no unused gene found in parent2, find any unused index
        if (!foundUnused) {
            for (size_t val = 0; val < size; ++val) {
                if (used.find(val) == used.end()) {
                    child[curr] = val;
                    used.insert(val);
                    foundUnused = true;
                    break;
                }
            }
        }
        
        // If still no unused gene found, something is wrong with our logic
        if (!foundUnused) {
            LOG_ERROR("OrderCrossover: Failed to find unused gene");
            return {};  // Return empty vector to signal failure
        }
        
        curr = (curr + 1) % size;
    }
    
    return child;
}

std::vector<std::shared_ptr<Individual>> PMXCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !parents[0] || !parents[1]) {
        LOG_WARNING("PMXCrossover: Invalid parents");
        return parents;
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    
    if (parent1Genes.empty() || parent2Genes.empty()) {
        LOG_WARNING("PMXCrossover: Empty parent genes");
        return parents;
    }
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);
    
    int maxAttempts = 3;
    while (offspring.empty() && maxAttempts > 0) {
        try {
            auto& rng = Random::instance();
            int point1 = rng.getRandomInt(0, static_cast<int>(parent1Genes.size() - 1));
            int point2 = rng.getRandomInt(0, static_cast<int>(parent1Genes.size() - 1));
            if (point1 > point2) std::swap(point1, point2);
            
            // Create offspring genes
            std::vector<int> child1Genes(parent1Genes.size());
            std::vector<int> child2Genes(parent2Genes.size());
            
            // Copy the mapping section
            for (int i = point1; i <= point2; i++) {
                child1Genes[i] = parent2Genes[i];
                child2Genes[i] = parent1Genes[i];
            }
            
            // Create mapping between values in the mapping section
            std::unordered_map<int, int> mapping1, mapping2;
            for (int i = point1; i <= point2; i++) {
                mapping1[parent1Genes[i]] = parent2Genes[i];
                mapping2[parent2Genes[i]] = parent1Genes[i];
            }
            
            // Fill remaining positions
            for (int i = 0; i < static_cast<int>(parent1Genes.size()); i++) {
                if (i >= point1 && i <= point2) continue;
                
                // For child1
                int value1 = parent1Genes[i];
                while (mapping1.find(value1) != mapping1.end()) {
                    value1 = mapping1[value1];
                }
                child1Genes[i] = value1;
                
                // For child2
                int value2 = parent2Genes[i];
                while (mapping2.find(value2) != mapping2.end()) {
                    value2 = mapping2[value2];
                }
                child2Genes[i] = value2;
            }
            
            // Create and validate offspring
            auto offspring1 = std::make_shared<Individual>(child1Genes);
            auto offspring2 = std::make_shared<Individual>(child2Genes);
            
            if (!offspring1->getGenes().empty() && representation->isValid(offspring1, instance)) {
                offspring.push_back(offspring1);
            }
            if (!offspring2->getGenes().empty() && representation->isValid(offspring2, instance)) {
                offspring.push_back(offspring2);
            }
            
        } catch (const std::exception& e) {
            LOG_WARNING("PMXCrossover: Exception during crossover: " + std::string(e.what()));
        }
        maxAttempts--;
    }
    
    // If we couldn't create any valid offspring, return copies of parents
    if (offspring.empty()) {
        LOG_WARNING("PMXCrossover: Failed to create valid offspring, returning parents");
        return parents;
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
            
            int overlapQuality = dna_utils::calculateEdgeWeight(
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
        // Log parent details
        LOG_DEBUG("Parent 1 genes size: " + std::to_string(parents[0]->getGenes().size()));
        LOG_DEBUG("Parent 2 genes size: " + std::to_string(parents[1]->getGenes().size()));
        
        // Log parent validity
        LOG_DEBUG("Parent 1 valid: " + std::to_string(representation->isValid(parents[0], instance)));
        LOG_DEBUG("Parent 2 valid: " + std::to_string(representation->isValid(parents[1], instance)));
        
        auto [genes1, genes2] = validateAndGetGenes(parents);
        if (genes1.empty() || genes2.empty()) {
            LOG_ERROR("Empty genes after validation");
            return parents;
        }
        
        // Find alignment segments in both parents
        auto segments1 = findAlignmentSegments(genes1, instance);
        auto segments2 = findAlignmentSegments(genes2, instance);
        
        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);
        
        // Create two offspring with different segment combinations
        for (int i = 0; i < 2; i++) {
            auto childGenes = mergeSegments(segments1, segments2, genes1.size());
            auto child = std::make_shared<Individual>(std::move(childGenes));
            if (!child->getGenes().empty() && representation->isValid(child, instance)) {
                offspring.push_back(child);
            } else {
                LOG_WARNING("Created invalid child in DNA alignment crossover");
            }
            
            // Swap segment orders for second child
            std::swap(segments1, segments2);
        }
        
        if (offspring.empty()) {
            LOG_WARNING("No valid offspring produced in DNA alignment crossover");
            return parents;
        }
        
        return offspring;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in DNA alignment crossover: " + std::string(e.what()));
        return parents;
    }
}

// Validate parents and get their genes
std::pair<std::vector<int>, std::vector<int>> validateAndGetGenes(
    const std::vector<std::shared_ptr<Individual>>& parents) {
    
    if (parents.size() < 2 || !parents[0] || !parents[1]) {
        LOG_ERROR("Invalid parents: size=" + std::to_string(parents.size()) + 
                 ", p1=" + (parents[0] ? "valid" : "null") + 
                 ", p2=" + (parents[1] ? "valid" : "null"));
        return {};
    }

    const auto& parent1 = parents[0];
    const auto& parent2 = parents[1];

    const auto& genes1 = parent1->getGenes();
    const auto& genes2 = parent2->getGenes();

    LOG_DEBUG("Validating genes - p1 size: " + std::to_string(genes1.size()) + 
              ", p2 size: " + std::to_string(genes2.size()));

    if (genes1.empty() || genes2.empty()) {
        LOG_ERROR("Empty genes in parents");
        return {};
    }

    if (genes1.size() != genes2.size()) {
        LOG_ERROR("Parents have different lengths: " + std::to_string(genes1.size()) + 
                 " vs " + std::to_string(genes2.size()));
        return {};
    }
    
    return {genes1, genes2};
} 