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
    if (!individual) return false;
    
    return validateGenes(individual->getGenes(), instance);
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
    if (!individual) {
        LOG_ERROR("Null individual provided to toDNA");
        return "";
    }
    
    const auto& genes = individual->getGenes();
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    
    if (genes.empty() || spectrum.empty()) {
        LOG_ERROR("Empty genes or spectrum in toDNA");
        return "";
    }
    
    if (k <= 1) {
        LOG_ERROR("Invalid k value: " + std::to_string(k));
        return "";
    }

    std::string dna;
    dna.reserve(instance.getN() * 2);  // Reserve extra space for safety
    
    // Find first valid k-mer to start with
    bool foundFirst = false;
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) {
            LOG_ERROR("Invalid gene index: " + std::to_string(genes[i]));
            continue;
        }
        
        try {
            const std::string& kmer = spectrum[genes[i]];
            if (kmer.length() >= static_cast<size_t>(k)) {
                dna = kmer;
                foundFirst = true;
                break;
            } else {
                LOG_ERROR("K-mer too short at index " + std::to_string(genes[i]) + 
                         ": " + std::to_string(kmer.length()) + " < " + std::to_string(k));
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Error accessing k-mer: " + std::string(e.what()));
            continue;
        }
    }
    
    if (!foundFirst) {
        LOG_ERROR("No valid starting k-mer found");
        return "";
    }
    
    // Add subsequent k-mers with proper overlap checking
    size_t maxAttempts = genes.size() * 2;  // Prevent infinite loops
    size_t attempts = 0;
    
    for (size_t i = 1; i < genes.size() && attempts < maxAttempts; ++i) {
        attempts++;
        
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) {
            LOG_ERROR("Invalid gene index: " + std::to_string(genes[i]));
            continue;
        }
        
        try {
            const std::string& kmer = spectrum[genes[i]];
            if (kmer.length() < static_cast<size_t>(k)) {
                LOG_ERROR("K-mer too short: " + kmer);
                continue;
            }
            
            // Verify overlap length is valid
            size_t overlapLen = k - 1;
            if (dna.length() < overlapLen || kmer.length() < overlapLen) {
                LOG_ERROR("Invalid overlap length - DNA: " + std::to_string(dna.length()) + 
                         ", k-mer: " + std::to_string(kmer.length()));
                continue;
            }
            
            // Verify overlap matches before adding
            std::string suffix = dna.substr(dna.length() - overlapLen);
            std::string prefix = kmer.substr(0, overlapLen);
            
            if (suffix.length() != overlapLen || prefix.length() != overlapLen) {
                LOG_ERROR("Invalid overlap segment lengths - suffix: " + std::to_string(suffix.length()) + 
                         ", prefix: " + std::to_string(prefix.length()));
                continue;
            }
            
            // Only add if overlap is good (within deltaK mismatches)
            int mismatches = 0;
            for (size_t j = 0; j < overlapLen; ++j) {
                if (suffix[j] != prefix[j]) mismatches++;
            }
            
            if (mismatches <= instance.getDeltaK()) {
                std::string toAdd = kmer.substr(overlapLen);
                if (!toAdd.empty()) {
                    dna += toAdd;
                } else {
                    LOG_ERROR("Empty k-mer segment to add");
                }
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Error processing k-mer overlap: " + std::string(e.what()));
            continue;
        }
    }
    
    if (dna.empty()) {
        LOG_ERROR("Failed to construct DNA sequence");
        return "";
    }
    
    if (dna.length() < static_cast<size_t>(k)) {
        LOG_ERROR("Constructed DNA sequence too short: " + std::to_string(dna.length()));
        return "";
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
        if (gene < 0 || gene >= static_cast<int>(spectrum.size())) {
            return false;
        }
        if (!used[gene] && spectrum[gene].length() >= static_cast<size_t>(instance.getK())) {
            used[gene] = true;
            usedValidKmers++;
        }
    }
    
    // Require at least 70% of valid k-mers to be used
    return usedValidKmers >= static_cast<int>(validKmers * 0.7);
} 