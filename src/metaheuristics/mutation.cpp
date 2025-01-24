//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/mutation.h"
#include "utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>

// =========== PointMutation ===========

void PointMutation::mutate(
    std::vector<std::shared_ptr<std::vector<int>>>& offspring,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (!offspring.empty() && !offspring[0]) {
        LOG_ERROR("Null individual passed to mutation");
        return;
    }

    LOG_DEBUG("Starting mutation on offspring of size " + std::to_string(offspring.size()));
    
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> distProb(0.0, 1.0);
    std::uniform_int_distribution<int> distPos(0, (int)offspring[0]->size() - 1);
    
    int mutations = 0;
    for (auto &ind : offspring) {
        if (!ind || ind->size() < 2) {
            continue;
        }
        double p = distProb(rng);
        if (p < mutationRate) {
            int posA = distPos(rng);
            int posB = distPos(rng);
            if (posA != posB) {
                std::swap((*ind)[posA], (*ind)[posB]);
                mutations++;
                
                DEBUG_LOG("Mutated position " + std::to_string(posA) + 
                         " from " + std::to_string((*ind)[posA]) + 
                         " to " + std::to_string((*ind)[posB]));
            }
        }
    }
    
    if (mutations > 0) {
        LOG_INFO("Applied " + std::to_string(mutations) + " mutations to offspring");
    } else {
        LOG_DEBUG("No mutations applied to offspring");
    }
}
