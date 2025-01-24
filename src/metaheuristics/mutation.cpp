//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "../include/metaheuristics/mutation_impl.h"
#include "../include/utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>

// =========== PointMutation ===========

void PointMutation::mutate(
    std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation
) {
    if (!solution || solution->empty()) {
        return;
    }
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    // Mutate each position with probability m_mutationRate
    for (size_t i = 0; i < solution->size(); i++) {
        if (dis(gen) < m_mutationRate) {
            // Flip bit or change value based on representation
            (*solution)[i] = !(*solution)[i];
        }
    }
}
