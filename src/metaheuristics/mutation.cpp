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

PointMutation::PointMutation(double mutationRate) : m_mutationRate(mutationRate) {
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        throw std::invalid_argument("Mutation rate must be between 0 and 1");
    }
}

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual) {
        LOG_ERROR("Null individual in mutation");
        return;
    }
    
    auto& genes = individual->getGenes();
    if (genes.size() < 2) {
        LOG_ERROR("Individual has too few genes for mutation");
        return;
    }
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    
    // For each position, attempt mutation with probability m_mutationRate
    for (size_t i = 0; i < genes.size(); ++i) {
        if (dist(gen) < m_mutationRate) {
            // Select random position to swap with
            std::uniform_int_distribution<size_t> posDist(0, genes.size() - 1);
            size_t j = posDist(gen);
            
            // Swap genes
            std::swap(genes[i], genes[j]);
        }
    }
    
    // Create new individual with mutated genes
    auto mutated = std::make_shared<Individual>(genes);
    
    // Only apply mutation if it results in a valid individual
    if (representation->isValid(mutated, instance)) {
        individual = mutated;
    }
}

// Implementation moved to mutation_impl.cpp
