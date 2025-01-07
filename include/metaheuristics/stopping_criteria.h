//
// Created by konrad_guest on 07/01/2025.
//

#ifndef STOPPING_CRITERIA_H
#define STOPPING_CRITERIA_H

#include <vector>
#include "metaheuristics/fitness.h"
#include "generator/dna_generator.h"
#include "metaheuristics/representation.h"

class IStopping {
public:
    virtual ~IStopping() = default;

    virtual bool stop(const std::vector<void*>         &population,
                      int                               generation,
                      const DNAInstance                &instance,
                      std::shared_ptr<IFitness>         fitness,
                      std::shared_ptr<IRepresentation>  representation) = 0;
};

/** Example: Stop by max generations */
class MaxGenerationsStopping : public IStopping {
public:
    explicit MaxGenerationsStopping(int maxGen) : m_maxGen(maxGen) {}

    bool stop(const std::vector<void*>         &population,
              int                               generation,
              const DNAInstance                &instance,
              std::shared_ptr<IFitness>         fitness,
              std::shared_ptr<IRepresentation>  representation) override
    {
        return generation >= m_maxGen;
    }

private:
    int m_maxGen;
};

// class NoImprovementStopping : public IStopping {
// public:
//     NoImprovementStopping(int maxNoImprove);
//     bool stop(const std::vector<std::vector<double>> &population,
//               int generation,
//               const IFitness &fitness) override;
// private:
//     int m_maxNoImprove;
//     double m_bestSoFar;
//     int m_noImproveCount;
// };
//
// class TimeLimitStopping : public IStopping {
// public:
//     TimeLimitStopping(double limitSec);
//     bool stop(const std::vector<std::vector<double>> &population,
//               int generation,
//               const IFitness &fitness) override;
// private:
//     double m_limitSec;
//     double m_startTime;
// };

#endif //STOPPING_CRITERIA_H
