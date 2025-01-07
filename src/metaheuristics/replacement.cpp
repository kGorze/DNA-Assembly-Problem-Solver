//
// Created by konrad_guest on 07/01/2025.
//

#include "metaheuristics/replacement.h"
#include <algorithm>

// =========== FullReplacement ===========
std::vector<void*> 
FullReplacement::replace(const std::vector<void*>        &oldPop,
                         const std::vector<void*>        &offspring,
                         const DNAInstance               &instance,
                         std::shared_ptr<IFitness>        fitness,
                         std::shared_ptr<IRepresentation> representation)
{
    // trivial: discard oldPop, return offspring
    return offspring;
}

// // =========== PartialReplacement ===========
// std::vector<std::vector<double>> 
// PartialReplacement::replace(const std::vector<std::vector<double>> &oldPop,
//                             const std::vector<std::vector<double>> &offspring,
//                             const IFitness &fitness)
// {
//     //TODO
//     // - zastąpić część populacji offspring
//     // ...
//     return {};
// }
//
// // =========== SteadyStateReplacement ===========
// std::vector<std::vector<double>> 
// SteadyStateReplacement::replace(const std::vector<std::vector<double>> &oldPop,
//                                 const std::vector<std::vector<double>> &offspring,
//                                 const IFitness &fitness)
// {
//     // ...
//     return {};
// }