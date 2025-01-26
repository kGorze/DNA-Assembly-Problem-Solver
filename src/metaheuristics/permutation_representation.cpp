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
    std::vector<int> genes;
    genes.resize(instance.getSpectrum().size());
    std::iota(genes.begin(), genes.end(), 0);  // Fill with 0, 1, 2, ..., n-1
    
    // Shuffle using random indices
    for (size_t i = genes.size() - 1; i > 0; --i) {
        int j = rng.getRandomInt(0, static_cast<int>(i));
        std::swap(genes[i], genes[j]);
    }
    
    individual.setGenes(genes);
    return true;
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
    if (!individual) return false;
    
    return validateGenes(individual->getGenes(), instance);
}

std::string PermutationRepresentation::toString(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance) const {
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
    
    // Add first k-mer completely
    if (!genes.empty() && genes[0] < static_cast<int>(spectrum.size())) {
        dna = spectrum[genes[0]];
    }
    
    // Add subsequent k-mers with proper overlap
    for (size_t i = 1; i < genes.size(); ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) continue;
        
        const std::string& kmer = spectrum[genes[i]];
        if (kmer.length() >= static_cast<size_t>(instance.getK())) {
            // Add only the non-overlapping part of each k-mer
            dna += kmer.substr(instance.getK() - 1);
        }
    }
    
    return dna;
}

bool PermutationRepresentation::validateGenes(
    const std::vector<int>& genes,
    const DNAInstance& instance) const {
    if (genes.empty()) return false;

    const size_t spectrumSize = instance.getSpectrum    ().size();
    std::vector<bool> used(spectrumSize, false);
    
    // Check bounds and mark used indices
    for (int gene : genes) {
        if (gene < 0 || gene >= static_cast<int>(spectrumSize)) {
            return false;
        }
        used[gene] = true;
    }
    
    // Check if all required indices are used (relaxed validation)
    int usedCount = std::count(used.begin(), used.end(), true);
    return usedCount >= static_cast<int>(spectrumSize * 0.9);  // Allow up to 10% missing indices
} 