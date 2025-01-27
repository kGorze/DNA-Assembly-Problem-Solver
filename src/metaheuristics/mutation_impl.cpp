#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/utils/logging.h"
#include <memory>
#include <random>
#include <stdexcept>
#include <sstream>

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual) {
        LOG_ERROR("Individual is null");
        return;
    }

    auto genes = individual->getGenes();
    if (genes.size() < 2) {
        LOG_ERROR("Individual has less than 2 genes, cannot perform mutation");
        return;
    }

    // Create a copy for mutation
    auto mutated = std::make_shared<Individual>(genes);
    auto& mutatedGenes = mutated->getGenes();
    
    // Calculate number of mutations based on size, ensuring minimum mutations
    int numMutations = std::max(m_minMutations, static_cast<int>(genes.size() * m_mutationRate));
    
    auto& rng = Random::instance();
    bool anyValidMutation = false;
    int consecutiveFailures = 0;
    
    // Try multiple mutations with backtracking
    for (int i = 0; i < numMutations && consecutiveFailures < 5; ++i) {
        // Store current state
        auto currentState = mutatedGenes;
        
        // Try up to 3 different positions for a successful mutation
        bool mutationSuccessful = false;
        for (int attempt = 0; attempt < 3 && !mutationSuccessful; ++attempt) {
            // Select two random positions for swapping
            int pos1 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
            int pos2;
            do {
                pos2 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
            } while (pos1 == pos2);
            
            // Perform swap
            std::swap(mutatedGenes[pos1], mutatedGenes[pos2]);
            
            // Try additional random swap to increase diversity
            if (rng.generateProbability() < 0.3) {  // 30% chance for additional swap
                int pos3 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
                int pos4;
                do {
                    pos4 = rng.getRandomInt(0, static_cast<int>(mutatedGenes.size() - 1));
                } while (pos3 == pos4);
                std::swap(mutatedGenes[pos3], mutatedGenes[pos4]);
            }
            
            // Check if this mutation made any improvement
            auto tempIndividual = std::make_shared<Individual>(mutatedGenes);
            if (representation->isValid(tempIndividual, instance)) {
                mutationSuccessful = true;
                anyValidMutation = true;
                consecutiveFailures = 0;
            } else {
                // Undo this mutation and try another position
                mutatedGenes = currentState;
            }
        }
        
        if (!mutationSuccessful) {
            consecutiveFailures++;
            mutatedGenes = currentState;  // Revert to last valid state
        }
    }
    
    // Only update if we made valid changes and the result is different
    if (anyValidMutation && mutatedGenes != genes) {
        individual = mutated;
        LOG_DEBUG("PointMutation: Successfully performed multiple mutations");
    } else {
        LOG_DEBUG("PointMutation: No valid mutations found or no changes made");
    }
} 