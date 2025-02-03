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
#include <climits>
#include <sstream>
#include <deque>

// Forward declaration of helper functions
std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size);
bool validateGeneSequence(const std::vector<int>& genes, const DNAInstance& instance);
int calculateKmerOverlapQuality(const std::vector<int>& genes, const DNAInstance& instance, size_t start, size_t length);

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
        
        if (size <= 1) {
            LOG_ERROR("Cannot generate random points: size must be greater than 1");
            return parents;
        }
        int point1 = rng.getRandomInt(0, static_cast<int>(size - 1));
        int point2 = rng.getRandomInt(0, static_cast<int>(size - 1));
        if (point1 < 0 || point2 < 0 || point1 >= static_cast<int>(size) || point2 >= static_cast<int>(size)) {
            LOG_ERROR("Generated random points out of bounds: " + std::to_string(point1) + ", " + std::to_string(point2));
            return parents;
        }
        if (point1 > point2) std::swap(point1, point2);
        
        int maxAttempts = 3;
        while (offspring.size() < 2 && maxAttempts > 0) {
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
        return parents;
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    const size_t spectrumSize = instance.getSpectrum().size();
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    
    // Validate gene lengths and values
    if (parent1Genes.empty() || parent2Genes.empty() || 
        parent1Genes.size() != spectrumSize || parent2Genes.size() != spectrumSize) {
        LOG_WARNING("EdgeRecombination: Invalid parent genes size");
        return parents;
    }
    
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(2);
    
    try {
        // Build edge table with overlap weights and local structure
        std::unordered_map<int, std::vector<std::pair<int, double>>> edgeTable;
        for (size_t i = 0; i < spectrumSize; i++) {
            edgeTable[i] = std::vector<std::pair<int, double>>();
        }
        
        // Add edges from both parents with overlap weights and position-based bonuses
        auto addEdges = [&](const std::vector<int>& genes, [[maybe_unused]] const char* parentName) {
            for (size_t i = 0; i < genes.size(); i++) {
                int current = genes[i];
                // Look at neighbors within a window of 5 positions
                for (int offset = -5; offset <= 5; offset++) {
                    if (offset == 0) continue;
                    size_t j = (i + genes.size() + offset) % genes.size();
                    int neighbor = genes[j];
                    
                    // Calculate base overlap quality
                    double overlapScore = dna_utils::calculateEdgeWeight(
                        spectrum[current], spectrum[neighbor], k);
                    
                    // Add position-based bonuses
                    if (std::abs(offset) == 1) {
                        overlapScore *= 2.0;  // Double weight for adjacent positions
                    } else if (std::abs(offset) <= 3) {
                        overlapScore *= 1.5;  // 50% bonus for near positions
                    }
                    
                    if (overlapScore > 0) {
                        edgeTable[current].push_back({neighbor, overlapScore});
                    }
                }
            }
        };
        
        addEdges(parent1Genes, "Parent1");
        addEdges(parent2Genes, "Parent2");
        
        // Create two offspring
        for (int childIdx = 0; childIdx < 2; childIdx++) {
            std::vector<int> childGenes;
            childGenes.reserve(spectrumSize);
            std::vector<bool> used(spectrumSize, false);
            
            // Start with a gene that has good overlaps from the respective parent
            const auto& startGenes = (childIdx == 0) ? parent1Genes : parent2Genes;
            int currentGene = startGenes[0];
            
            // Track the last few genes for local structure preservation
            std::deque<int> recentGenes;
            const size_t MEMORY_SIZE = 5;
            
            childGenes.push_back(currentGene);
            used[currentGene] = true;
            recentGenes.push_back(currentGene);
            
            while (childGenes.size() < spectrumSize) {
                int nextGene = -1;
                double bestQuality = -1;
                
                // Get candidates with their qualities
                std::vector<std::pair<int, double>> candidates;
                
                // First try: look at neighbors with direct overlaps
                if (!edgeTable[currentGene].empty()) {
                    for (const auto& [neighbor, quality] : edgeTable[currentGene]) {
                        if (!used[neighbor]) {
                            double adjustedQuality = quality;
                            
                            // Check overlap with recent genes
                            if (!recentGenes.empty()) {
                                for (size_t i = 0; i < recentGenes.size(); i++) {
                                    int recentGene = recentGenes[i];
                                    double recentOverlap = dna_utils::calculateEdgeWeight(
                                        spectrum[neighbor], spectrum[recentGene], k);
                                    adjustedQuality += recentOverlap * (1.0 - i * 0.2);  // Decay factor
                                }
                            }
                            
                            candidates.push_back({neighbor, adjustedQuality});
                        }
                    }
                }
                
                // Sort candidates by quality
                std::sort(candidates.begin(), candidates.end(),
                         [](const auto& a, const auto& b) { return a.second > b.second; });
                
                // Try top candidates
                const int MAX_CANDIDATES = 5;
                for (size_t i = 0; i < std::min(candidates.size(), size_t(MAX_CANDIDATES)); i++) {
                    int candidate = candidates[i].first;
                    // Verify k-mer overlap quality
                    if (!childGenes.empty()) {
                        int overlap = dna_utils::calculateEdgeWeight(
                            spectrum[currentGene], spectrum[candidate], k);
                        if (overlap > 0) {
                            nextGene = candidate;
                            bestQuality = candidates[i].second;
                            break;
                        }
                    } else {
                        nextGene = candidate;
                        bestQuality = candidates[i].second;
                        break;
                    }
                }
                
                // If no good candidate found, try to repair
                if (nextGene == -1) {
                    // Look for any unused gene with acceptable overlap
                    for (size_t i = 0; i < spectrumSize; i++) {
                        if (!used[i]) {
                            int overlap = dna_utils::calculateEdgeWeight(
                                spectrum[currentGene], spectrum[i], k);
                            if (overlap > bestQuality) {
                                bestQuality = overlap;
                                nextGene = i;
                            }
                        }
                    }
                }
                
                // If still no gene found, take first unused
                if (nextGene == -1) {
                    for (size_t i = 0; i < spectrumSize; i++) {
                        if (!used[i]) {
                            nextGene = i;
                            break;
                        }
                    }
                }
                
                childGenes.push_back(nextGene);
                used[nextGene] = true;
                currentGene = nextGene;
                
                // Update recent genes
                recentGenes.push_back(nextGene);
                if (recentGenes.size() > MEMORY_SIZE) {
                    recentGenes.pop_front();
                }
            }
            
            // Attempt local repair of poor overlaps
            bool improved;
            int repairAttempts = 0;
            const int MAX_REPAIR_ATTEMPTS = 5;
            
            do {
                improved = false;
                for (size_t i = 0; i < childGenes.size() - 1; i++) {
                    int curr = childGenes[i];
                    int next = childGenes[i + 1];
                    int overlap = dna_utils::calculateEdgeWeight(spectrum[curr], spectrum[next], k);
                    
                    if (overlap == 0) {  // Poor overlap found
                        // Try to find a better gene to insert between them
                        int bestInsert = -1;
                        int bestCombinedOverlap = -1;
                        
                        for (size_t j = 0; j < spectrumSize; j++) {
                            if (j != static_cast<size_t>(curr) && j != static_cast<size_t>(next)) {
                                int overlapPrev = dna_utils::calculateEdgeWeight(spectrum[curr], spectrum[j], k);
                                int overlapNext = dna_utils::calculateEdgeWeight(spectrum[j], spectrum[next], k);
                                int combinedOverlap = overlapPrev + overlapNext;
                                
                                if (combinedOverlap > bestCombinedOverlap) {
                                    bestCombinedOverlap = combinedOverlap;
                                    bestInsert = j;
                                }
                            }
                        }
                        
                        if (bestInsert != -1 && bestCombinedOverlap > 0) {
                            // Insert the new gene and remove it from its old position
                            auto oldPos = std::find(childGenes.begin(), childGenes.end(), bestInsert);
                            if (oldPos != childGenes.end()) {
                                childGenes.erase(oldPos);
                            }
                            childGenes.insert(childGenes.begin() + i + 1, bestInsert);
                            improved = true;
                        }
                    }
                }
                repairAttempts++;
            } while (improved && repairAttempts < MAX_REPAIR_ATTEMPTS);
            
            // Create and validate offspring
            auto child = std::make_shared<Individual>(childGenes);
            if (validateGeneSequence(childGenes, instance) && representation->isValid(child, instance)) {
                offspring.push_back(child);
                LOG_DEBUG("EdgeRecombination: Created valid offspring");
            } else {
                LOG_WARNING("EdgeRecombination: Child validation failed");
                // Try to use the better parent as fallback
                if (parents[childIdx] && representation->isValid(parents[childIdx], instance)) {
                    offspring.push_back(parents[childIdx]);
                    LOG_DEBUG("EdgeRecombination: Using parent as fallback");
                }
            }
        }
        
    } catch (const std::exception& e) {
        LOG_WARNING("EdgeRecombination: Exception during crossover: " + std::string(e.what()));
    }
    
    if (offspring.empty()) {
        LOG_WARNING("EdgeRecombination: Failed to create valid offspring");
        return parents;
    }
    
    return offspring;
}

std::vector<std::shared_ptr<Individual>> CycleCrossover::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    // Input validation
    if (parents.size() < 2 || !parents[0] || !parents[1] || !representation) {
        LOG_WARNING("CycleCrossover: Invalid input parameters");
        return {};
    }
    
    const auto& parent1Genes = parents[0]->getGenes();
    const auto& parent2Genes = parents[1]->getGenes();
    
    // Validate parent genes
    if (parent1Genes.empty() || parent2Genes.empty() || 
        parent1Genes.size() != parent2Genes.size() ||
        parent1Genes.size() != instance.getSpectrum().size()) {
        LOG_WARNING("CycleCrossover: Invalid parent genes");
        return {};
    }
    
    // Validate gene values
    const size_t size = parent1Genes.size();
    for (size_t i = 0; i < size; i++) {
        if (parent1Genes[i] < 0 || static_cast<size_t>(parent1Genes[i]) >= size ||
            parent2Genes[i] < 0 || static_cast<size_t>(parent2Genes[i]) >= size) {
            LOG_WARNING("CycleCrossover: Invalid gene values in parents");
            return {};
        }
    }
    
    // Create offspring with same size as parents
    std::vector<int> offspring1(size, -1);
    std::vector<int> offspring2(size, -1);
    std::vector<bool> used1(size, false);
    std::vector<bool> used2(size, false);
    
    // Find cycles and create offspring
    for (size_t start = 0; start < size; start++) {
        if (offspring1[start] != -1) continue;  // Position already filled
        
        // Find cycle starting at this position
        size_t pos = start;
        bool useParent1 = true;  // Alternate for each cycle
        
        while (offspring1[pos] == -1) {
            // Copy genes from appropriate parent based on cycle
            if (useParent1) {
                offspring1[pos] = parent1Genes[pos];
                offspring2[pos] = parent2Genes[pos];
            } else {
                offspring1[pos] = parent2Genes[pos];
                offspring2[pos] = parent1Genes[pos];
            }
            
            used1[offspring1[pos]] = true;
            used2[offspring2[pos]] = true;
            
            // Find next position in cycle
            int value = parent2Genes[pos];
            pos = std::find(parent1Genes.begin(), parent1Genes.end(), value) - parent1Genes.begin();
            
            if (pos >= size || offspring1[pos] != -1) break;  // End of cycle or already filled
        }
        
        useParent1 = !useParent1;  // Switch for next cycle
    }
    
    // Validate offspring
    auto validateOffspring = [&](const std::vector<int>& genes) -> bool {
        std::vector<bool> seen(size, false);
        for (int gene : genes) {
            if (gene < 0 || static_cast<size_t>(gene) >= size || seen[gene]) {
                return false;
            }
            seen[gene] = true;
        }
        return std::find(seen.begin(), seen.end(), false) == seen.end();
    };
    
    std::vector<std::shared_ptr<Individual>> result;
    
    if (validateOffspring(offspring1)) {
        auto child1 = std::make_shared<Individual>(offspring1);
        if (representation->isValid(child1, instance)) {
            result.push_back(child1);
        }
    }
    
    if (validateOffspring(offspring2)) {
        auto child2 = std::make_shared<Individual>(offspring2);
        if (representation->isValid(child2, instance)) {
            result.push_back(child2);
        }
    }
    
    // If no valid offspring were created, return empty vector
    if (result.empty()) {
        LOG_WARNING("CycleCrossover: Failed to create valid offspring");
        return {};
    }
    
    return result;
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
            
            // Create offspring genes
            std::vector<int> offspring1Genes(size, -1);
            std::vector<int> offspring2Genes(size, -1);
            
            // Step 1: Copy the mapping section
            for (int i = point1; i <= point2; i++) {
                offspring1Genes[i] = parent2Genes[i];
                offspring2Genes[i] = parent1Genes[i];
            }
            
            // Step 2: Create mapping between values in the mapping section
            std::unordered_map<int, int> mapping1, mapping2;
            for (int i = point1; i <= point2; i++) {
                mapping1[parent1Genes[i]] = parent2Genes[i];
                mapping2[parent2Genes[i]] = parent1Genes[i];
            }
            
            // Step 3: Fill remaining positions
            for (int i = 0; i < static_cast<int>(size); i++) {
                if (i >= point1 && i <= point2) continue;
                
                // For offspring1
                int value1 = parent1Genes[i];
                while (true) {
                    auto it = mapping1.find(value1);
                    if (it == mapping1.end()) break;
                    value1 = it->second;
                }
                offspring1Genes[i] = value1;
                
                // For offspring2
                int value2 = parent2Genes[i];
                while (true) {
                    auto it = mapping2.find(value2);
                    if (it == mapping2.end()) break;
                    value2 = it->second;
                }
                offspring2Genes[i] = value2;
            }
            
            // Validate offspring before creating individuals
            auto validateOffspring = [&](const std::vector<int>& genes) -> bool {
                if (genes.size() != size) return false;
                std::vector<bool> used(size, false);
                for (int gene : genes) {
                    if (gene < 0 || static_cast<size_t>(gene) >= size || used[gene]) {
                        return false;
                    }
                    used[gene] = true;
                }
                return true;
            };
            
            if (validateOffspring(offspring1Genes)) {
                auto child1 = std::make_shared<Individual>(offspring1Genes);
                if (representation->isValid(child1, instance)) {
                    offspring.push_back(child1);
                }
            }
            
            if (validateOffspring(offspring2Genes)) {
                auto child2 = std::make_shared<Individual>(offspring2Genes);
                if (representation->isValid(child2, instance)) {
                    offspring.push_back(child2);
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

std::vector<int> performOrderCrossover(const std::vector<int>& parent1, const std::vector<int>& parent2, size_t size) {
    // Input validation
    if (parent1.empty() || parent2.empty() || size == 0 || 
        size > parent1.size() || size > parent2.size()) {
        LOG_ERROR("OrderCrossover: Invalid input parameters");
        return {};
    }

    // Validate parent genes are within valid range
    for (int gene : parent1) {
        if (gene < 0 || static_cast<size_t>(gene) >= size) {
            LOG_ERROR("OrderCrossover: Parent 1 contains invalid gene value: " + std::to_string(gene));
            return {};
        }
    }
    for (int gene : parent2) {
        if (gene < 0 || static_cast<size_t>(gene) >= size) {
            LOG_ERROR("OrderCrossover: Parent 2 contains invalid gene value: " + std::to_string(gene));
            return {};
        }
    }
    
    auto& rng = Random::instance();
    
    // Select two random crossover points
    int point1 = rng.getRandomInt(0, static_cast<int>(size - 1));
    int point2 = rng.getRandomInt(0, static_cast<int>(size - 1));
    if (point1 > point2) std::swap(point1, point2);
    
    // Initialize child with -1 to mark unfilled positions
    std::vector<int> child(size, -1);
    std::vector<bool> used(size, false);
    
    // Step 1: Copy segment from parent1
    for (int i = point1; i <= point2; i++) {
        child[i] = parent1[i];
        used[parent1[i]] = true;
    }
    
    // Step 2: Fill remaining positions with genes from parent2 in order
    int curr_pos = (point2 + 1) % size;  // Start filling from position after segment
    int p2_pos = (point2 + 1) % size;    // Start taking genes from parent2 after segment
    
    while (curr_pos != point1) {
        // Find next unused gene from parent2
        while (p2_pos != point2 + 1) {
            int candidate = parent2[p2_pos];
            if (!used[candidate]) {
                child[curr_pos] = candidate;
                used[candidate] = true;
                break;
            }
            p2_pos = (p2_pos + 1) % size;
        }
        
        // If we haven't filled this position yet, find first unused value
        if (child[curr_pos] == -1) {
            for (size_t val = 0; val < size; val++) {
                if (!used[val]) {
                    child[curr_pos] = static_cast<int>(val);
                    used[val] = true;
                    break;
                }
            }
        }
        
        curr_pos = (curr_pos + 1) % size;
    }
    
    // Validate the final child
    for (size_t i = 0; i < size; i++) {
        if (child[i] == -1 || child[i] < 0 || static_cast<size_t>(child[i]) >= size) {
            LOG_ERROR("OrderCrossover: Invalid child generated");
            return {};
        }
    }
    
    // Verify it's a valid permutation
    std::vector<bool> verify_used(size, false);
    for (int gene : child) {
        if (verify_used[gene]) {
            LOG_ERROR("OrderCrossover: Generated child contains duplicate genes");
            return {};
        }
        verify_used[gene] = true;
    }
    
    return child;
}

bool validateGeneSequence(const std::vector<int>& genes, const DNAInstance& instance) {
    if (genes.empty()) {
        LOG_ERROR("Empty gene sequence");
        return false;
    }

    // Check for size match with spectrum
    if (genes.size() != instance.getSpectrum().size()) {
        LOG_ERROR("Gene sequence size mismatch: " + std::to_string(genes.size()) + 
                 " vs expected " + std::to_string(instance.getSpectrum().size()));
        return false;
    }

    // Check for valid gene indices and duplicates
    std::vector<bool> used(instance.getSpectrum().size(), false);
    for (int gene : genes) {
        if (gene < 0 || static_cast<size_t>(gene) >= instance.getSpectrum().size()) {
            LOG_ERROR("Invalid gene index: " + std::to_string(gene));
            return false;
        }
        if (used[gene]) {
            LOG_ERROR("Duplicate gene found: " + std::to_string(gene));
            return false;
        }
        used[gene] = true;
    }

    // Check k-mer overlaps
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    for (size_t i = 0; i < genes.size() - 1; i++) {
        int curr = genes[i];
        int next = genes[i + 1];
        int overlap = dna_utils::calculateEdgeWeight(spectrum[curr], spectrum[next], k);
        if (overlap <= 0) {
            LOG_ERROR("Invalid k-mer overlap at position " + std::to_string(i) + 
                     ": " + spectrum[curr] + " -> " + spectrum[next]);
            return false;
        }
    }

    return true;
}

int calculateKmerOverlapQuality(const std::vector<int>& genes, 
                              const DNAInstance& instance,
                              size_t start,
                              size_t length) {
    if (start + length > genes.size()) return 0;
    
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    int totalOverlap = 0;
    
    for (size_t i = start; i < start + length - 1; i++) {
        int curr = genes[i];
        int next = genes[i + 1];
        if (curr >= 0 && next >= 0 && 
            static_cast<size_t>(curr) < spectrum.size() && 
            static_cast<size_t>(next) < spectrum.size()) {
            int overlap = dna_utils::calculateEdgeWeight(spectrum[curr], spectrum[next], k);
            totalOverlap += overlap;
        }
    }
    
    return totalOverlap;
} 