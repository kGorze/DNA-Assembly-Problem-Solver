#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/metaheuristics/edge_table.h"
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

namespace {

// Utility functions for handling duplicate genes
[[maybe_unused]]
bool hasDuplicateGenes(const std::vector<int>& genes) {
    std::unordered_set<int> geneSet;
    for (int gene : genes) {
        if (gene >= 0 && !geneSet.insert(gene).second) {
            return true;
        }
    }
    return false;
}

[[maybe_unused]]
void repairDuplicateGenes(std::vector<int>& genes, const DNAInstance& instance) {
    const size_t spectrumSize = instance.getSpectrum().size();
    std::vector<bool> used(spectrumSize, false);
    std::vector<size_t> duplicatePositions;
    
    // Mark used genes and find duplicates
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] >= 0 && static_cast<size_t>(genes[i]) < spectrumSize) {
            if (used[genes[i]]) {
                duplicatePositions.push_back(i);
            } else {
                used[genes[i]] = true;
            }
        } else {
            duplicatePositions.push_back(i);
        }
    }
    
    // Replace duplicates with unused genes
    for (size_t pos : duplicatePositions) {
        for (size_t i = 0; i < spectrumSize; ++i) {
            if (!used[i]) {
                genes[pos] = static_cast<int>(i);
                used[i] = true;
                break;
            }
        }
    }
}

} // anonymous namespace

std::vector<std::shared_ptr<Individual>> OrderCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    try {
        auto [genes1, genes2] = validateAndGetGenes(parents);
        if (genes1.empty() || genes2.empty()) {
            LOG_WARNING("Invalid parent genes in order crossover");
            return parents;
        }
        
        std::vector<std::shared_ptr<Individual>> offspring;
        offspring.reserve(2);
        
        auto& rng = Random::instance();
        const size_t size = std::min(genes1.size(), genes2.size());
        
        int maxAttempts = 3;
        while (offspring.size() < 2 && maxAttempts > 0) {
            int start = rng.getRandomInt(0, static_cast<int>(size - 1));
            int end = rng.getRandomInt(0, static_cast<int>(size - 1));
            if (start > end) std::swap(start, end);
            
            // Create offspring with duplicate detection and repair
            auto offspring1Genes = performOrderCrossover(genes1, genes2, size);
            if (!offspring1Genes.empty()) {
                if (auto child = createOffspring(offspring1Genes, representation, instance)) {
                    offspring.push_back(child);
                }
            }
            
            auto offspring2Genes = performOrderCrossover(genes2, genes1, size);
            if (!offspring2Genes.empty()) {
                if (auto child = createOffspring(offspring2Genes, representation, instance)) {
                    offspring.push_back(child);
                }
            }
            
            maxAttempts--;
        }
        
        if (offspring.empty()) {
            LOG_WARNING("Failed to create valid offspring in order crossover, validating parents before return");
            std::vector<std::shared_ptr<Individual>> validatedParents;
            for (const auto& parent : parents) {
                if (parent && representation->isValid(parent, instance)) {
                    validatedParents.push_back(parent);
                }
            }
            if (!validatedParents.empty()) {
                return validatedParents;
            }
            LOG_ERROR("OrderCrossover: No valid parents to return as fallback");
            return {};
        }
        
        return offspring;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in order crossover: " + std::string(e.what()));
        return parents;
    }
}

std::vector<std::shared_ptr<Individual>> EdgeRecombinationCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (parents.size() < 2 || !parents[0] || !parents[1]) {
        LOG_WARNING("EdgeRecombination: Invalid parents");
        return {};
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    const size_t spectrumSize = instance.getSpectrum().size();
    
    // Validate gene lengths and values
    if (parent1Genes.empty() || parent2Genes.empty()) {
        LOG_WARNING("EdgeRecombination: Empty parent genes");
        return {};
    }
    
    // Enforce same length equal to spectrum size
    if (parent1Genes.size() != spectrumSize || parent2Genes.size() != spectrumSize) {
        LOG_WARNING("EdgeRecombination: Parents must have length equal to spectrum size (" + 
                   std::to_string(spectrumSize) + ")");
        return {};
    }
    
    // Validate gene values in parents
    for (const auto& genes : {parent1Genes, parent2Genes}) {
        for (int gene : genes) {
            if (gene < 0 || static_cast<size_t>(gene) >= spectrumSize) {
                LOG_WARNING("EdgeRecombination: Invalid gene value in parent");
                return {};
            }
        }
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
                childGenes.reserve(spectrumSize);
                
                // Start with a random gene from first parent
                int startPos = rng.getRandomInt(0, static_cast<int>(spectrumSize - 1));
                if (startPos < 0 || static_cast<size_t>(startPos) >= parent1Genes.size()) {
                    LOG_ERROR("EdgeRecombination: Invalid start position generated");
                    continue;
                }
                
                int current = parent1Genes[startPos];
                if (current < 0 || static_cast<size_t>(current) >= spectrumSize) {
                    LOG_ERROR("EdgeRecombination: Invalid start gene selected");
                    continue;
                }
                
                childGenes.push_back(current);
                std::unordered_set<int> usedGenes{current};
                
                // Build the rest of the path prioritizing high-quality overlaps
                bool validPath = true;
                while (childGenes.size() < spectrumSize && validPath) {
                    const auto neighbors = edgeTable.getNeighbors(current);
                    bool foundNext = false;
                    
                    // First try: Find unused neighbor with best overlap quality
                    for (const auto& neighbor : neighbors) {
                        if (neighbor.node < 0 || static_cast<size_t>(neighbor.node) >= spectrumSize) {
                            LOG_ERROR("EdgeRecombination: Invalid neighbor node");
                            validPath = false;
                            break;
                        }
                        
                        if (usedGenes.find(neighbor.node) == usedGenes.end()) {
                            childGenes.push_back(neighbor.node);
                            usedGenes.insert(neighbor.node);
                            current = neighbor.node;
                            foundNext = true;
                            break;
                        }
                    }
                    
                    if (!foundNext && validPath) {
                        // Second try: Use any unused gene from parents
                        std::vector<int> unusedGenes;
                        for (int gene : parent1Genes) {
                            if (gene >= 0 && static_cast<size_t>(gene) < spectrumSize && 
                                usedGenes.find(gene) == usedGenes.end()) {
                                unusedGenes.push_back(gene);
                            }
                        }
                        
                        if (!unusedGenes.empty()) {
                            int idx = rng.getRandomInt(0, static_cast<int>(unusedGenes.size() - 1));
                            if (idx >= 0 && idx < static_cast<int>(unusedGenes.size())) {
                                int nextGene = unusedGenes[idx];
                                if (nextGene >= 0 && static_cast<size_t>(nextGene) < spectrumSize) {
                                    childGenes.push_back(nextGene);
                                    usedGenes.insert(nextGene);
                                    current = nextGene;
                                    foundNext = true;
                                }
                            }
                        }
                        
                        if (!foundNext) {
                            validPath = false;
                        }
                    }
                }
                
                // Validate final child genes
                bool validChild = true;
                if (validPath && childGenes.size() == spectrumSize) {
                    for (int gene : childGenes) {
                        if (gene < 0 || static_cast<size_t>(gene) >= spectrumSize) {
                            validChild = false;
                            break;
                        }
                    }
                    
                    if (validChild) {
                        auto child = std::make_shared<Individual>(childGenes);
                        if (representation->isValid(child, instance)) {
                            offspring.push_back(child);
                            LOG_DEBUG("EdgeRecombination: Created valid offspring");
                        }
                    }
                }
            }
            maxAttempts--;
            
        } catch (const std::exception& e) {
            LOG_WARNING("EdgeRecombination: Exception during crossover: " + std::string(e.what()));
            maxAttempts--;
        }
    }
    
    if (offspring.empty()) {
        LOG_WARNING("EdgeRecombination: Failed to create valid offspring, validating parents before return");
        std::vector<std::shared_ptr<Individual>> validatedParents;
        for (const auto& parent : parents) {
            if (parent && representation->isValid(parent, instance)) {
                validatedParents.push_back(parent);
            }
        }
        if (!validatedParents.empty()) {
            return validatedParents;
        }
        LOG_ERROR("EdgeRecombination: No valid parents to return as fallback");
        return {};
    }
    
    return offspring;
}

std::vector<std::shared_ptr<Individual>> CycleCrossover::crossover(
    [[maybe_unused]] const std::vector<std::shared_ptr<Individual>>& parents,
    [[maybe_unused]] const DNAInstance& instance,
    [[maybe_unused]] std::shared_ptr<IRepresentation> representation) {
    
    // TODO: Implement cycle crossover
    LOG_WARNING("CycleCrossover not implemented yet");
    return parents;  // Return parents as fallback
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
    
    // Use the smaller parent size to avoid buffer overflows
    size_t size = std::min(parent1Genes.size(), parent2Genes.size());
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);
    
    int maxAttempts = 3;
    while (offspring.empty() && maxAttempts > 0) {
        try {
            auto& rng = Random::instance();
            int point1 = rng.getRandomInt(0, static_cast<int>(size - 1));
            int point2 = rng.getRandomInt(0, static_cast<int>(size - 1));
            if (point1 > point2) std::swap(point1, point2);
            
            // Create offspring genes with proper size
            std::vector<int> child1Genes(size, -1);  // Initialize with -1 to detect unset positions
            std::vector<int> child2Genes(size, -1);
            
            // Copy the mapping section
            for (int i = point1; i <= point2; i++) {
                if (i < static_cast<int>(parent2Genes.size())) {
                    child1Genes[i] = parent2Genes[i];
                }
                if (i < static_cast<int>(parent1Genes.size())) {
                    child2Genes[i] = parent1Genes[i];
                }
            }
            
            // Create mapping between values in the mapping section
            std::unordered_map<int, int> mapping1, mapping2;
            for (int i = point1; i <= point2; i++) {
                if (i < static_cast<int>(std::min(parent1Genes.size(), parent2Genes.size()))) {
                    mapping1[parent1Genes[i]] = parent2Genes[i];
                    mapping2[parent2Genes[i]] = parent1Genes[i];
                }
            }
            
            // Fill remaining positions
            for (int i = 0; i < static_cast<int>(size); i++) {
                if (i >= point1 && i <= point2) continue;
                
                // For child1
                if (i < static_cast<int>(parent1Genes.size())) {
                    int value1 = parent1Genes[i];
                    int iterations = 0;
                    std::unordered_set<int> seen;
                    while (mapping1.find(value1) != mapping1.end() && iterations < 100) {
                        if (seen.find(value1) != seen.end()) {
                            // Cycle detected, use first unused value
                            for (size_t val = 0; val < instance.getSpectrum().size(); val++) {
                                if (mapping1.find(val) == mapping1.end()) {
                                    value1 = val;
                                    break;
                                }
                            }
                            break;
                        }
                        seen.insert(value1);
                        value1 = mapping1[value1];
                        iterations++;
                    }
                    child1Genes[i] = value1;
                }
                
                // For child2
                if (i < static_cast<int>(parent2Genes.size())) {
                    int value2 = parent2Genes[i];
                    int iterations = 0;
                    std::unordered_set<int> seen;
                    while (mapping2.find(value2) != mapping2.end() && iterations < 100) {
                        if (seen.find(value2) != seen.end()) {
                            // Cycle detected, use first unused value
                            for (size_t val = 0; val < instance.getSpectrum().size(); val++) {
                                if (mapping2.find(val) == mapping2.end()) {
                                    value2 = val;
                                    break;
                                }
                            }
                            break;
                        }
                        seen.insert(value2);
                        value2 = mapping2[value2];
                        iterations++;
                    }
                    child2Genes[i] = value2;
                }
            }
            
            // Fill any remaining -1 positions with unused valid indices
            std::vector<bool> used1(instance.getSpectrum().size(), false);
            std::vector<bool> used2(instance.getSpectrum().size(), false);
            
            // Mark used values
            for (int gene : child1Genes) {
                if (gene >= 0 && static_cast<size_t>(gene) < instance.getSpectrum().size()) {
                    used1[gene] = true;
                }
            }
            for (int gene : child2Genes) {
                if (gene >= 0 && static_cast<size_t>(gene) < instance.getSpectrum().size()) {
                    used2[gene] = true;
                }
            }
            
            // Fill remaining positions
            for (size_t i = 0; i < size; i++) {
                // Fill child1
                if (child1Genes[i] == -1) {
                    for (size_t val = 0; val < instance.getSpectrum().size(); val++) {
                        if (!used1[val]) {
                            child1Genes[i] = static_cast<int>(val);
                            used1[val] = true;
                            break;
                        }
                    }
                }
                
                // Fill child2
                if (child2Genes[i] == -1) {
                    for (size_t val = 0; val < instance.getSpectrum().size(); val++) {
                        if (!used2[val]) {
                            child2Genes[i] = static_cast<int>(val);
                            used2[val] = true;
                            break;
                        }
                    }
                }
            }
            
            // Validate final offspring
            bool valid1 = true, valid2 = true;
            for (int gene : child1Genes) {
                if (gene < 0 || static_cast<size_t>(gene) >= instance.getSpectrum().size()) {
                    valid1 = false;
                    break;
                }
            }
            for (int gene : child2Genes) {
                if (gene < 0 || static_cast<size_t>(gene) >= instance.getSpectrum().size()) {
                    valid2 = false;
                    break;
                }
            }
            
            // Create and validate offspring
            if (valid1) {
                auto offspring1 = std::make_shared<Individual>(child1Genes);
                if (representation->isValid(offspring1, instance)) {
                    offspring.push_back(offspring1);
                }
            }
            if (valid2) {
                auto offspring2 = std::make_shared<Individual>(child2Genes);
                if (representation->isValid(offspring2, instance)) {
                    offspring.push_back(offspring2);
                }
            }
            
        } catch (const std::exception& e) {
            LOG_WARNING("PMXCrossover: Exception during crossover: " + std::string(e.what()));
        }
        maxAttempts--;
    }
    
    if (offspring.empty()) {
        LOG_WARNING("PMXCrossover: Failed to create valid offspring, validating parents before return");
        std::vector<std::shared_ptr<Individual>> validatedParents;
        for (const auto& parent : parents) {
            if (parent && representation->isValid(parent, instance)) {
                validatedParents.push_back(parent);
            }
        }
        if (!validatedParents.empty()) {
            return validatedParents;
        }
        LOG_ERROR("PMXCrossover: No valid parents to return as fallback");
        return {};
    }
    
    return offspring;
}

std::vector<DNAAlignmentCrossover::AlignmentSegment> DNAAlignmentCrossover::findAlignmentSegments(
    const std::vector<int>& genes,
    const DNAInstance& instance,
    int minSegmentLength) const {
    
    std::vector<AlignmentSegment> segments;
    
    // Find all possible segments of length >= minLength
    for (size_t start = 0; start < genes.size(); start++) {
        for (size_t len = minSegmentLength; start + len <= genes.size(); len++) {
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

std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size) {
    // Use the shared Random instance instead of creating a new one
    auto& rng = Random::instance();
    
    // Validate size
    if (size == 0 || size > parent1.size() || size > parent2.size()) {
        LOG_ERROR("OrderCrossover: Invalid size parameter");
        return {};
    }
    
    // Select two random crossover points with proper bounds checking
    int point1 = rng.getRandomInt(0, static_cast<int>(size - 1));
    int point2 = rng.getRandomInt(0, static_cast<int>(size - 1));
    if (point1 > point2) std::swap(point1, point2);
    
    // Validate crossover points
    if (point1 < 0 || point2 < 0 || static_cast<size_t>(point1) >= size || static_cast<size_t>(point2) >= size) {
        LOG_ERROR("OrderCrossover: Invalid crossover points generated");
        return {};
    }
    
    // Create child with same size as parents
    std::vector<int> child(size);
    std::unordered_set<int> used;
    
    // Copy segment between crossover points from parent1 with validation
    for (int i = point1; i <= point2; ++i) {
        if (static_cast<size_t>(i) >= parent1.size()) {
            LOG_ERROR("OrderCrossover: Index out of bounds while copying segment");
            return {};
        }
        child[i] = parent1[i];
        used.insert(parent1[i]);
    }
    
    // Fill remaining positions with genes from parent2 in order
    int curr = (point2 + 1) % size;
    int p2pos = (point2 + 1) % size;
    
    while (curr != point1) {
        // Try to find an unused gene from parent2
        bool foundUnused = false;
        for (size_t i = 0; i < size; ++i) {
            int checkPos = (p2pos + static_cast<int>(i)) % size;
            if (checkPos < 0 || static_cast<size_t>(checkPos) >= parent2.size()) {
                LOG_ERROR("OrderCrossover: Invalid position while checking parent2");
                return {};
            }
            
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
                if (used.find(static_cast<int>(val)) == used.end()) {
                    child[curr] = static_cast<int>(val);
                    used.insert(static_cast<int>(val));
                    foundUnused = true;
                    break;
                }
            }
        }
        
        // If still no unused gene found, something is wrong
        if (!foundUnused) {
            LOG_ERROR("OrderCrossover: Failed to find unused gene");
            return {};
        }
        
        curr = (curr + 1) % size;
    }
    
    // Validate final child
    for (size_t i = 0; i < child.size(); ++i) {
        if (child[i] < 0 || static_cast<size_t>(child[i]) >= size) {
            LOG_ERROR("OrderCrossover: Invalid gene value in final child");
            return {};
        }
    }
    
    return child;
} 