//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <sstream>

// =========== PointMutation ===========

PointMutation::PointMutation(double mutationRate) {
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        throw std::invalid_argument("Mutation rate must be between 0 and 1");
    }
    m_mutationRate = mutationRate;
}

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual) {
        LOG_ERROR("Individual is null");
        return;
    }

    if (individual->getSize() < 2) {
        LOG_ERROR("Individual has less than 2 genes, cannot perform mutation");
        return;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> posDis(0, individual->getSize() - 1);

    if (dis(gen) < m_mutationRate) {
        int i = posDis(gen);
        int j;
        do {
            j = posDis(gen);
        } while (i == j);

        // Create a copy of the individual for mutation
        auto mutated = individual->clone();
        auto& genes = mutated->getGenes();  // Get non-const reference
        std::swap(genes[i], genes[j]);

        // Only apply mutation if it results in a valid individual
        if (representation->isValid(mutated, instance)) {
            individual = mutated;
        }
    }
}

// Implementation moved to mutation_impl.cpp
