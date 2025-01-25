#pragma once

#include "interfaces/i_fitness.h"
#include "../dna/dna_instance.h"
#include "individual.h"
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <mutex>

class SimpleFitness : public IFitness {
private:
    static constexpr size_t maxCacheSize = 10000;
    std::unordered_map<std::string, double> cache;
    std::mutex cacheMutex;
    
    double calculateConnectivityScore(const std::vector<int>& genes,
                                    const DNAInstance& instance) const;
    double calculateSpectrumCoverageScore(const std::vector<int>& genes,
                                        const DNAInstance& instance) const;
    double calculateLengthPenalty(const std::vector<int>& genes,
                                const DNAInstance& instance) const;
    void cleanupCache();
    
public:
    SimpleFitness() = default;
    ~SimpleFitness() override = default;
    
    double evaluate(const std::shared_ptr<Individual>& individual,
                   const DNAInstance& instance) override;
}; 