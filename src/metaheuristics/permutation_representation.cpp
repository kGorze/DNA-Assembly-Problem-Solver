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
    
    // First, identify valid k-mers and their overlap properties
    std::vector<int> validKmers;
    std::vector<std::vector<int>> goodOverlaps(spectrum.size());  // Store indices of k-mers that overlap well
    
    // Collect valid k-mers
    for (size_t i = 0; i < spectrum.size(); ++i) {
        if (spectrum[i].length() >= static_cast<size_t>(instance.getK())) {
            validKmers.push_back(i);
        }
    }
    
    // Build overlap graph
    for (size_t i = 0; i < validKmers.size(); ++i) {
        int idx1 = validKmers[i];
        const std::string& kmer1 = spectrum[idx1];
        
        for (size_t j = 0; j < validKmers.size(); ++j) {
            if (i == j) continue;
            
            int idx2 = validKmers[j];
            const std::string& kmer2 = spectrum[idx2];
            
            // Check suffix-prefix overlap
            std::string suffix = kmer1.substr(kmer1.length() - (instance.getK() - 1));
            std::string prefix = kmer2.substr(0, instance.getK() - 1);
            
            int mismatches = 0;
            for (size_t k = 0; k < static_cast<size_t>(instance.getK() - 1); ++k) {
                if (suffix[k] != prefix[k]) mismatches++;
            }
            
            // Store good overlaps (allow one more mismatch than deltaK for flexibility)
            if (mismatches <= instance.getDeltaK() + 1) {
                goodOverlaps[idx1].push_back(idx2);
            }
        }
    }
    
    // Start with a random valid k-mer that has good overlaps
    std::vector<int> selectedGenes;
    std::vector<bool> used(spectrum.size(), false);
    
    // Find starting k-mer with good connectivity
    std::vector<int> goodStarters;
    for (int idx : validKmers) {
        if (goodOverlaps[idx].size() >= 2) {  // Has at least 2 good overlaps
            goodStarters.push_back(idx);
        }
    }
    
    if (goodStarters.empty()) {
        // Fall back to any valid k-mer if no good starters found
        if (validKmers.empty()) return false;
        selectedGenes.push_back(validKmers[rng.getRandomInt(0, static_cast<int>(validKmers.size()) - 1)]);
    } else {
        selectedGenes.push_back(goodStarters[rng.getRandomInt(0, static_cast<int>(goodStarters.size()) - 1)]);
    }
    used[selectedGenes[0]] = true;
    
    // Grow the sequence by preferring good overlaps
    while (selectedGenes.size() < validKmers.size() * 0.7) {  // Try to use at least 70% of valid k-mers
        int lastIdx = selectedGenes.back();
        const auto& possibleNext = goodOverlaps[lastIdx];
        
        std::vector<int> candidates;
        // Prefer unused k-mers with good overlaps
        for (int nextIdx : possibleNext) {
            if (!used[nextIdx]) {
                candidates.push_back(nextIdx);
            }
        }
        
        if (candidates.empty()) {
            // If no good overlaps available, try any unused valid k-mer
            for (int idx : validKmers) {
                if (!used[idx]) {
                    candidates.push_back(idx);
                }
            }
            
            if (candidates.empty()) break;  // No more unused k-mers
        }
        
        // Select next k-mer
        int nextIdx = candidates[rng.getRandomInt(0, static_cast<int>(candidates.size()) - 1)];
        selectedGenes.push_back(nextIdx);
        used[nextIdx] = true;
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
    if (!individual) return "";
    
    const auto& genes = individual->getGenes();
    const auto& spectrum = instance.getSpectrum();
    if (genes.empty() || spectrum.empty()) return "";

    std::string dna;
    dna.reserve(instance.getN() * 2);  // Reserve extra space for safety
    
    // Find first valid k-mer to start with
    bool foundFirst = false;
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) continue;
        
        const std::string& kmer = spectrum[genes[i]];
        if (kmer.length() >= static_cast<size_t>(instance.getK())) {
            dna = kmer;
            foundFirst = true;
            break;
        }
    }
    
    if (!foundFirst) return "";
    
    // Add subsequent k-mers with proper overlap checking
    for (size_t i = 1; i < genes.size(); ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) continue;
        
        const std::string& kmer = spectrum[genes[i]];
        if (kmer.length() < static_cast<size_t>(instance.getK())) continue;
        
        // Verify overlap length is valid
        size_t overlapLen = instance.getK() - 1;
        if (dna.length() < overlapLen || kmer.length() < overlapLen) continue;
        
        // Verify overlap matches before adding
        std::string suffix = dna.substr(dna.length() - overlapLen);
        std::string prefix = kmer.substr(0, overlapLen);
        
        // Only add if overlap is good (within deltaK mismatches)
        int mismatches = 0;
        for (size_t j = 0; j < overlapLen; ++j) {
            if (suffix[j] != prefix[j]) mismatches++;
        }
        
        if (mismatches <= instance.getDeltaK()) {
            dna += kmer.substr(overlapLen);
        }
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