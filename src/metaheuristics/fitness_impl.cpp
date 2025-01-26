#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/utils/logging.h"

double SimpleFitness::calculateFitness(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) const {
    
    if (!solution || !representation) {
        LOG_ERROR("Invalid solution or representation in SimpleFitness");
        return 0.0;
    }

    // Calculate connectivity score
    double connectivityScore = calculateConnectivityScore(solution, instance);
    
    // Get DNA string representation
    std::string dnaStr = representation->toString(solution, instance);
    std::vector<char> dna(dnaStr.begin(), dnaStr.end());
    
    // Calculate spectrum coverage score
    double coverageScore = calculateSpectrumCoverageScore(dna, instance);
    
    // Calculate length penalty
    double lengthPenalty = calculateLengthPenalty(dna.size(), instance.getN());
    
    // Combine scores with weights
    return 0.4 * connectivityScore + 0.4 * coverageScore + 0.2 * lengthPenalty;
}

double SimpleFitness::calculateConnectivityScore(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    
    // Simple implementation - check if adjacent genes form valid connections
    const auto& genes = solution->getGenes();
    if (genes.empty()) return 0.0;
    
    int validConnections = 0;
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (std::abs(genes[i] - genes[i + 1]) <= 1) {
            validConnections++;
        }
    }
    
    return static_cast<double>(validConnections) / (genes.size() - 1);
}

double SimpleFitness::calculateSpectrumCoverageScore(
    const std::vector<char>& dna,
    const DNAInstance& instance) const {
    
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty() || dna.empty()) return 0.0;
    
    int k = instance.getK();
    int covered = 0;
    
    for (const auto& oligo : spectrum) {
        bool found = false;
        for (size_t i = 0; i <= dna.size() - k; ++i) {
            std::string subseq(dna.begin() + i, dna.begin() + i + k);
            if (subseq == oligo) {
                found = true;
                break;
            }
        }
        if (found) covered++;
    }
    
    return static_cast<double>(covered) / spectrum.size();
}

double SimpleFitness::calculateLengthPenalty(
    int actualLength,
    int targetLength) const {
    
    if (actualLength == targetLength) return 1.0;
    
    double diff = std::abs(actualLength - targetLength);
    return std::max(0.0, 1.0 - diff / targetLength);
} 