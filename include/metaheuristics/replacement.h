//
// Created by konrad_guest on 07/01/2025.
//

#ifndef REPLACEMENT_H
#define REPLACEMENT_H

#include <vector>
#include "metaheuristics/fitness.h"

/**
 * Strategia zastępowania pokoleń:
 * - całościowe zastępowanie
 * - częściowe zastępowanie
 * - steady-state
 */
class IReplacement {
public:
    virtual ~IReplacement() = default;
    virtual std::vector<void*>
    replace(const std::vector<void*>        &oldPop,
                         const std::vector<void*>        &offspring,
                         const DNAInstance               &instance,
                         std::shared_ptr<IFitness>        fitness,
                         std::shared_ptr<IRepresentation> representation) = 0;
};

class FullReplacement : public IReplacement {
public:
    virtual std::vector<void*>
    replace(const std::vector<void*>        &oldPop,
                         const std::vector<void*>        &offspring,
                         const DNAInstance               &instance,
                         std::shared_ptr<IFitness>        fitness,
                         std::shared_ptr<IRepresentation> representation) override;
};

// class PartialReplacement : public IReplacement {
// public:
//     std::vector<std::vector<double>> 
//     replace(const std::vector<std::vector<double>> &oldPop,
//             const std::vector<std::vector<double>> &offspring,
//             const IFitness &fitness) override;
// };
//
// class SteadyStateReplacement : public IReplacement {
// public:
//     std::vector<std::vector<double>> 
//     replace(const std::vector<std::vector<double>> &oldPop,
//             const std::vector<std::vector<double>> &offspring,
//             const IFitness &fitness) override;
// };

#endif //REPLACEMENT_H
