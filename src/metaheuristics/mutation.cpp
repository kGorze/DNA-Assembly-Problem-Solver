//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "../include/metaheuristics/mutation_impl.h"
#include "../include/utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>

// =========== PointMutation ===========

void PointMutation::mutate(
    std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) {
    if (!solution || solution->empty()) {
        LOG_WARNING("Attempted to mutate null or empty solution");
        return;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to mutation operator");
        return;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> posDis(0, solution->size() - 1);
    std::uniform_int_distribution<> valDis(0, instance.getSpectrum().size() - 1);

    auto originalSolution = *solution;
    bool mutated = false;

    for (size_t i = 0; i < solution->size(); ++i) {
        if (dis(gen) < m_mutationRate) {
            int oldValue = (*solution)[i];
            int newValue;
            do {
                newValue = valDis(gen);
            } while (newValue == oldValue);
            
            (*solution)[i] = newValue;
            mutated = true;

            std::stringstream ss;
            ss << "Mutated position " << i << " from " << oldValue << " to " << newValue;
            LOG_DEBUG(ss.str());
        }
    }

    if (!mutated) {
        // Force at least one mutation if none occurred
        int pos = posDis(gen);
        int oldValue = (*solution)[pos];
        int newValue;
        do {
            newValue = valDis(gen);
        } while (newValue == oldValue);
        
        (*solution)[pos] = newValue;
        
        std::stringstream ss;
        ss << "Forced mutation at position " << pos << " from " << oldValue << " to " << newValue;
        LOG_DEBUG(ss.str());
    }

    // Validate the mutated solution
    auto mutatedSolution = std::make_shared<std::vector<int>>(*solution);
    if (!representation->isValid(mutatedSolution, instance)) {
        LOG_WARNING("Mutation produced invalid solution - rolling back changes");
        *solution = originalSolution;
    }
}
