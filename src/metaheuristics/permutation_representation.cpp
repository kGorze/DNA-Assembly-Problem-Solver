#include "../../include/metaheuristics/representation.h"
#include "../../include/utils/random.h"
#include <algorithm>
#include <numeric>

std::vector<std::shared_ptr<Individual>> PermutationRepresentation::initializePopulation(
    int populationSize, const DNAInstance& instance) {
    if (populationSize <= 0) {
        throw std::invalid_argument("Population size must be positive");
    }
    
    std::vector<std::shared_ptr<Individual>> population;
    population.reserve(populationSize);
    
    for (int i = 0; i < populationSize; ++i) {
        auto individual = std::make_shared<Individual>();
        if (initializeIndividual(*individual, instance)) {
            population.push_back(individual);
        }
    }
    
    return population;
}

bool PermutationRepresentation::initializeIndividual(Individual& individual, const DNAInstance& instance) {
    auto& rng = Random::instance();
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    
    if (k <= 1 || spectrum.empty()) {
        LOG_ERROR("Invalid k value or empty spectrum");
        return false;
    }
    
    // First, identify valid k-mers and their overlap properties
    std::vector<int> validKmers;
    std::vector<std::vector<int>> goodOverlaps(spectrum.size());  // Store indices of k-mers that overlap well
    
    // Collect valid k-mers with length checks
    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (spectrum[i].length() >= static_cast<size_t>(k)) {
            validKmers.push_back(i);
        } else {
            LOG_ERROR("K-mer at index " + std::to_string(i) + " is too short: " + 
                      std::to_string(spectrum[i].length()) + " < " + std::to_string(k));
        }
    }
    
    if (validKmers.empty()) {
        LOG_ERROR("No valid k-mers found in spectrum");
        return false;
    }
    
    // Build overlap graph with bounds checking
    for (size_t i = 0; i < validKmers.size(); ++i) {
        int idx1 = validKmers[i];
        if (idx1 < 0 || idx1 >= static_cast<int>(spectrum.size())) {
            LOG_ERROR("Invalid k-mer index: " + std::to_string(idx1));
            continue;
        }
        
        const std::string& kmer1 = spectrum[idx1];
        if (kmer1.length() < static_cast<size_t>(k - 1)) {
            LOG_ERROR("K-mer too short for overlap: " + std::to_string(kmer1.length()));
            continue;
        }
        
        for (size_t j = 0; j < validKmers.size(); ++j) {
            if (i == j) continue;
            
            int idx2 = validKmers[j];
            if (idx2 < 0 || idx2 >= static_cast<int>(spectrum.size())) {
                LOG_ERROR("Invalid k-mer index: " + std::to_string(idx2));
                continue;
            }
            
            const std::string& kmer2 = spectrum[idx2];
            if (kmer2.length() < static_cast<size_t>(k - 1)) {
                LOG_ERROR("K-mer too short for overlap: " + std::to_string(kmer2.length()));
                continue;
            }
            
            try {
                // Check suffix-prefix overlap with bounds checking
                std::string suffix = kmer1.substr(kmer1.length() - (k - 1));
                std::string prefix = kmer2.substr(0, k - 1);
                
                if (suffix.length() != static_cast<size_t>(k - 1) || prefix.length() != static_cast<size_t>(k - 1)) {
                    LOG_ERROR("Invalid overlap lengths: suffix=" + std::to_string(suffix.length()) + 
                             ", prefix=" + std::to_string(prefix.length()));
                    continue;
                }
                
                int mismatches = 0;
                for (size_t pos = 0; pos < static_cast<size_t>(k - 1); ++pos) {
                    if (suffix[pos] != prefix[pos]) mismatches++;
                }
                
                // Store good overlaps (allow one more mismatch than deltaK for flexibility)
                if (mismatches <= instance.getDeltaK() + 1) {
                    goodOverlaps[idx1].push_back(idx2);
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Error processing k-mer overlap: " + std::string(e.what()));
                continue;
            }
        }
    }
    
    // Start with a random valid k-mer that has good overlaps
    std::vector<int> selectedGenes;
    std::vector<bool> used(spectrum.size(), false);
    
    // Find starting k-mer with good connectivity
    std::vector<int> goodStarters;
    for (int idx : validKmers) {
        if (idx >= 0 && idx < static_cast<int>(goodOverlaps.size()) && goodOverlaps[idx].size() >= 2) {
            goodStarters.push_back(idx);
        }
    }
    
    if (goodStarters.empty()) {
        // Fall back to any valid k-mer if no good starters found
        if (validKmers.empty()) {
            LOG_ERROR("No valid k-mers available for initialization");
            return false;
        }
        int randomIdx = rng.getRandomInt(0, static_cast<int>(validKmers.size()) - 1);
        if (randomIdx >= 0 && randomIdx < static_cast<int>(validKmers.size())) {
            selectedGenes.push_back(validKmers[randomIdx]);
        } else {
            LOG_ERROR("Invalid random index generated: " + std::to_string(randomIdx));
            return false;
        }
    } else {
        int randomIdx = rng.getRandomInt(0, static_cast<int>(goodStarters.size()) - 1);
        if (randomIdx >= 0 && randomIdx < static_cast<int>(goodStarters.size())) {
            selectedGenes.push_back(goodStarters[randomIdx]);
        } else {
            LOG_ERROR("Invalid random index generated: " + std::to_string(randomIdx));
            return false;
        }
    }
    
    if (selectedGenes.empty()) {
        LOG_ERROR("Failed to select initial k-mer");
        return false;
    }
    
    used[selectedGenes[0]] = true;
    
    // Grow the sequence by preferring good overlaps
    size_t targetSize = static_cast<size_t>(validKmers.size() * 0.7);  // Try to use at least 70% of valid k-mers
    size_t maxAttempts = validKmers.size() * 2;  // Prevent infinite loops
    size_t attempts = 0;
    
    while (selectedGenes.size() < targetSize && attempts < maxAttempts) {
        attempts++;
        int lastIdx = selectedGenes.back();
        
        if (lastIdx < 0 || lastIdx >= static_cast<int>(goodOverlaps.size())) {
            LOG_ERROR("Invalid last index: " + std::to_string(lastIdx));
            break;
        }
        
        const auto& possibleNext = goodOverlaps[lastIdx];
        std::vector<int> candidates;
        
        // Prefer unused k-mers with good overlaps
        for (int nextIdx : possibleNext) {
            if (nextIdx >= 0 && nextIdx < static_cast<int>(spectrum.size()) && !used[nextIdx]) {
                candidates.push_back(nextIdx);
            }
        }
        
        if (candidates.empty()) {
            // If no good overlaps available, try any unused valid k-mer
            for (int idx : validKmers) {
                if (idx >= 0 && idx < static_cast<int>(spectrum.size()) && !used[idx]) {
                    candidates.push_back(idx);
                }
            }
            
            if (candidates.empty()) break;  // No more unused k-mers
        }
        
        // Select next k-mer
        int randomIdx = rng.getRandomInt(0, static_cast<int>(candidates.size()) - 1);
        if (randomIdx >= 0 && randomIdx < static_cast<int>(candidates.size())) {
            int nextIdx = candidates[randomIdx];
            if (nextIdx >= 0 && nextIdx < static_cast<int>(spectrum.size())) {
                selectedGenes.push_back(nextIdx);
                used[nextIdx] = true;
            } else {
                LOG_ERROR("Invalid next index selected: " + std::to_string(nextIdx));
                break;
            }
        } else {
            LOG_ERROR("Invalid random index generated: " + std::to_string(randomIdx));
            break;
        }
    }
    
    if (selectedGenes.size() < 2) {
        LOG_ERROR("Failed to generate enough genes: " + std::to_string(selectedGenes.size()));
        return false;
    }
    
    individual.setGenes(selectedGenes);
    return validateGenes(selectedGenes, instance);
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    
    if (!individual) {
        LOG_WARNING("Null individual in isValid check");
        return false;
    }

    const auto& genes = individual->getGenes();
    if (genes.empty()) {
        LOG_WARNING("Empty genes in isValid check");
        return false;
    }

    const auto& spectrum = instance.getSpectrum();
    
    // Basic validation - just check if indices are within bounds
    for (int gene : genes) {
        if (gene < 0 || gene >= static_cast<int>(spectrum.size())) {
            return false;
        }
    }
    
    return true;  // Accept any solution with valid indices
}

std::string PermutationRepresentation::toString(
    const std::shared_ptr<Individual>& individual,
    [[maybe_unused]] const DNAInstance& instance) const {
    if (!individual) return "";
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return "";

    std::string result;
    result.reserve(genes.size() * 3);  // Estimate 3 chars per number
    
    for (size_t i = 0; i < genes.size(); ++i) {
        if (i > 0) result += " ";
        result += std::to_string(genes[i]);
    }
    
    return result;
}

std::string PermutationRepresentation::toDNA(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    if (!individual) return "";
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) return "";
    
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    const int targetLength = instance.getN();
    std::string dna = spectrum[genes[0]];  // Start with first k-mer
    
    // Track used k-mers to avoid duplicates
    std::vector<bool> used(spectrum.size(), false);
    used[genes[0]] = true;
    
    // For each remaining gene
    for (size_t i = 1; i < genes.size(); i++) {
        if (static_cast<size_t>(genes[i]) >= spectrum.size()) continue;  // Skip invalid indices
        
        // Skip if already used and repetitions not allowed
        if (!instance.isRepAllowed() && used[genes[i]]) continue;
        
        const std::string& current = spectrum[genes[i]];
        
        // Try different overlap lengths
        int bestOverlap = -1;
        int minMismatches = k + 1;  // More than maximum possible
        
        // Try overlaps from k-1 down to 1
        for (int overlap = k - 1; overlap > 0; overlap--) {
            if (static_cast<int>(dna.length()) < overlap) continue;
            
            // Count mismatches in overlap region
            int mismatches = 0;
            for (int j = 0; j < overlap; j++) {
                if (dna[dna.length() - overlap + j] != current[j]) {
                    mismatches++;
                }
            }
            
            // If this is the best overlap so far, store it
            if (mismatches < minMismatches) {
                minMismatches = mismatches;
                bestOverlap = overlap;
            }
            
            // If we found a perfect match, no need to check smaller overlaps
            if (mismatches == 0) break;
        }
        
        // If we found any valid overlap, append with that overlap
        if (bestOverlap > 0) {
            // Even if there are mismatches, we'll still append but track the quality
            dna += current.substr(bestOverlap);
            used[genes[i]] = true;
        } else {
            // No overlap found - append with a single base overlap
            // This ensures we keep building the solution even with poor overlaps
            dna += current.substr(1);
            used[genes[i]] = true;
        }
        
        // If we've reached or exceeded the target length, stop
        if (static_cast<int>(dna.length()) >= targetLength) {
            break;
        }
    }
    
    // If the DNA is too short, extend it by appending remaining k-mers
    while (static_cast<int>(dna.length()) < targetLength && genes.size() > 0) {
        for (size_t i = 0; i < genes.size(); i++) {
            if (static_cast<size_t>(genes[i]) < spectrum.size() && (!used[genes[i]] || instance.isRepAllowed())) {
                dna += spectrum[genes[i]].substr(1);  // Append with minimal overlap
                used[genes[i]] = true;
                if (static_cast<int>(dna.length()) >= targetLength) break;
            }
        }
        // If we can't extend further, break to avoid infinite loop
        if (static_cast<int>(dna.length()) < targetLength && std::all_of(used.begin(), used.end(),
            [&](bool b) { return b || !instance.isRepAllowed(); })) {
            break;
        }
    }
    
    // Trim or pad the DNA to exactly match the target length
    if (static_cast<int>(dna.length()) > targetLength) {
        dna = dna.substr(0, targetLength);
    } else while (static_cast<int>(dna.length()) < targetLength) {
        dna += 'N';  // Pad with 'N' for missing bases
    }
    
    return dna;
}

bool PermutationRepresentation::validateGenes(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    if (genes.empty()) return false;

    const auto& spectrum = instance.getSpectrum();
    std::vector<bool> used(spectrum.size(), false);
    
    // Count valid k-mers
    int validKmers = 0;
    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (spectrum[i].length() >= static_cast<size_t>(instance.getK())) {
            validKmers++;
        }
    }
    
    // Check bounds and mark used indices
    int usedValidKmers = 0;
    for (int gene : genes) {
        if (gene >= 0 && gene < static_cast<int>(spectrum.size()) && 
            spectrum[gene].length() >= static_cast<size_t>(instance.getK())) {
            if (!used[gene]) {
                used[gene] = true;
                usedValidKmers++;
            }
        }
    }
    
    // More lenient validation - require at least 50% of valid k-mers to be used
    return usedValidKmers >= static_cast<int>(validKmers * 0.5);
} 