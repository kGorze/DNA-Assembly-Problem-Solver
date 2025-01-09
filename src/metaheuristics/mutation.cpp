//
// Created by konrad_guest on 07/01/2025.
//
#include "metaheuristics/mutation.h"

#include <random>
#include <algorithm>

// =========== PointMutation ===========
void PointMutation::mutate(std::vector<void*> &offspring,
                          const DNAInstance &instance,
                          std::shared_ptr<IRepresentation> representation)
{
    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> distProb(0.0, 1.0);

    for (auto &ind : offspring) {
        if (distProb(rng) < mutationRate) {
            auto permPtr = static_cast<std::vector<int>*>(ind);
            if (!permPtr || permPtr->size() < 2) continue;

            std::uniform_int_distribution<int> distPos(0, (int)permPtr->size()-1);
            int posA = distPos(rng);
            int posB = distPos(rng);
            if (posA != posB) {
                std::swap((*permPtr)[posA], (*permPtr)[posB]);
            }
        }
    }
}


// // =========== GaussianMutation ===========
// void GaussianMutation::mutate(std::vector<std::vector<double>> &offspring) {
//     // ...
// }
//
// // =========== InversionMutation ===========
// void InversionMutation::mutate(std::vector<std::vector<double>> &offspring) {
//     // ...
// }