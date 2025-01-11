//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/mutation.h"
#include <random>
#include <algorithm>
#include <iostream>

// =========== PointMutation ===========

void PointMutation::mutate(
    std::vector<std::shared_ptr<std::vector<int>>>& offspring,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> distProb(0.0, 1.0);

    for (auto &ind : offspring) {
        if (!ind || ind->size() < 2) {
            continue;
        }
        double p = distProb(rng);
        if (p < mutationRate) {
            std::uniform_int_distribution<int> distPos(0, (int)ind->size() - 1);
            int posA = distPos(rng);
            int posB = distPos(rng);
            if (posA != posB) {
                std::swap((*ind)[posA], (*ind)[posB]);
            }
        }
    }
}
