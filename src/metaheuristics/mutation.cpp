//
// Created by konrad_guest on 07/01/2025.
//
#include "metaheuristics/mutation.h"

#include <random>
#include <algorithm>

// =========== PointMutation ===========
void PointMutation::mutate(std::vector<void*>        &offspring,
                           const DNAInstance         &instance,
                           std::shared_ptr<IRepresentation> representation)
{
    // do nothing for now
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