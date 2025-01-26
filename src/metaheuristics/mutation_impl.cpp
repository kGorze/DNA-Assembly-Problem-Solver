#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/utils/logging.h"
#include <memory>
#include <random>
#include <stdexcept>
#include <sstream>

PointMutation::PointMutation(double mutationRate) {
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        throw std::invalid_argument("Mutation rate must be between 0 and 1");
    }
    m_mutationRate = mutationRate;
    LOG_INFO("PointMutation initialized with rate: " + std::to_string(mutationRate));
}

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual || !representation) {
        LOG_WARNING("Null individual or representation in mutation");
        return;
    }

    auto& genes = individual->getGenes();
    if (genes.empty()) {
        LOG_WARNING("Empty genes in mutation");
        return;
    }

    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_WARNING("Empty spectrum in mutation");
        return;
    }

    // Perform a few swaps without strict validation
    int numSwaps = std::max(1, static_cast<int>(genes.size() * 0.1)); // Swap up to 10% of genes
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> dis(0, genes.size() - 1);
    
    for (int i = 0; i < numSwaps; ++i) {
        size_t pos1 = dis(gen);
        size_t pos2 = dis(gen);
        
        // Basic bounds checking
        if (pos1 < genes.size() && pos2 < genes.size() && 
            genes[pos1] >= 0 && genes[pos1] < static_cast<int>(spectrum.size()) &&
            genes[pos2] >= 0 && genes[pos2] < static_cast<int>(spectrum.size())) {
            
            // Perform swap without validation
            std::swap(genes[pos1], genes[pos2]);
            LOG_DEBUG("Mutation: Swapped positions " + std::to_string(pos1) + 
                     " and " + std::to_string(pos2));
        }
    }
    
    // No validation or reverting - let fitness function handle penalties
} 