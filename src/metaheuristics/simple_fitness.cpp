#include "metaheuristics/simple_fitness.h"
#include "utils/logger.h"
#include <cmath>
#include <limits>

double SimpleFitness::evaluate(
    const std::shared_ptr<Individual>& individual,
    const DNAInstance& instance)
{
    if (!individual) {
        LOG_ERROR("Cannot evaluate null individual");
        return -std::numeric_limits<double>::infinity();
    }
    
    const auto& genes = individual->getGenes();
    if (genes.empty()) {
        LOG_ERROR("Cannot evaluate empty individual");
        return -std::numeric_limits<double>::infinity();
    }
    
    try {
        // Calculate hash for caching
        std::string key = individual->toString();
        
        // Check cache first
        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            auto it = cache.find(key);
            if (it != cache.end()) {
                return it->second;
            }
        }
        
        // Calculate fitness components
        double connectivityScore = calculateConnectivityScore(genes, instance);
        double spectrumCoverageScore = calculateSpectrumCoverageScore(genes, instance);
        double lengthPenalty = calculateLengthPenalty(genes, instance);
        
        // Combine scores with weights
        double fitness = (0.4 * connectivityScore + 
                        0.4 * spectrumCoverageScore - 
                        0.2 * lengthPenalty);
        
        // Cache result with thread safety
        {
            std::lock_guard<std::mutex> lock(cacheMutex);
            if (cache.size() >= maxCacheSize) {
                cleanupCache();
            }
            cache[key] = fitness;
        }
        
        return fitness;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Fitness evaluation failed: " + std::string(e.what()));
        return -std::numeric_limits<double>::infinity();
    }
}

double SimpleFitness::calculateConnectivityScore(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    if (genes.empty()) {
        return 0.0;
    }
    
    const auto& spectrum = instance.getSpectrum();
    size_t validConnections = 0;
    size_t totalPossibleConnections = genes.size() - 1;
    
    for (size_t i = 0; i < genes.size() - 1; ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size()) ||
            genes[i + 1] < 0 || genes[i + 1] >= static_cast<int>(spectrum.size())) {
            continue;
        }
        
        const std::string& kmer1 = spectrum[genes[i]];
        const std::string& kmer2 = spectrum[genes[i + 1]];
        
        if (kmer1.length() >= instance.getK() - 1 && 
            kmer2.length() >= instance.getK() - 1) {
            std::string suffix = kmer1.substr(kmer1.length() - (instance.getK() - 1));
            std::string prefix = kmer2.substr(0, instance.getK() - 1);
            
            if (suffix == prefix) {
                ++validConnections;
            }
        }
    }
    
    return totalPossibleConnections > 0 ? 
           static_cast<double>(validConnections) / totalPossibleConnections : 0.0;
}

double SimpleFitness::calculateSpectrumCoverageScore(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    if (genes.empty()) {
        return 0.0;
    }
    
    const auto& spectrum = instance.getSpectrum();
    std::vector<bool> covered(spectrum.size(), false);
    size_t coveredCount = 0;
    
    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i] < 0 || genes[i] >= static_cast<int>(spectrum.size())) {
            continue;
        }
        if (!covered[genes[i]]) {
            covered[genes[i]] = true;
            ++coveredCount;
        }
    }
    
    return spectrum.size() > 0 ? 
           static_cast<double>(coveredCount) / spectrum.size() : 0.0;
}

double SimpleFitness::calculateLengthPenalty(
    const std::vector<int>& genes,
    const DNAInstance& instance) const
{
    double expectedLength = static_cast<double>(instance.getN());
    double actualLength = static_cast<double>(genes.size());
    double lengthDiff = std::abs(actualLength - expectedLength);
    
    // Normalize penalty to [0,1] range using exponential decay
    return std::exp(-lengthDiff / expectedLength);
}

void SimpleFitness::cleanupCache() {
    if (cache.size() <= maxCacheSize / 2) return;
    
    std::vector<std::pair<std::string, double>> entries(
        cache.begin(), cache.end());
    
    // Keep only the most recent half
    size_t keepCount = maxCacheSize / 2;
    cache.clear();
    for (size_t i = entries.size() - keepCount; i < entries.size(); ++i) {
        cache.insert(entries[i]);
    }
} 