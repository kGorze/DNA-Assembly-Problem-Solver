//
// Created by konrad_guest on 08/01/2025.
// SMART

#pragma once

#include "../interfaces/i_benchmark.h"
#include "../interfaces/i_crossover.h"
#include "../interfaces/i_representation.h"
#include "../interfaces/i_fitness.h"
#include "../dna/dna_instance.h"
#include "../metaheuristics/crossover_impl.h"
#include "../metaheuristics/representation.h"
#include "../utils/logging.h"
#include "../utils/timer.h"
#include <memory>
#include <vector>
#include <string>

class CrossoverBenchmark : public IBenchmark {
public:
    CrossoverBenchmark(std::shared_ptr<ICrossover> crossover,
                      std::shared_ptr<IRepresentation> representation,
                      std::shared_ptr<IFitness> fitness);

    void runBenchmark(const DNAInstance& instance) override;
    std::string getName() const override;

private:
    std::shared_ptr<ICrossover> m_crossover;
    std::shared_ptr<IRepresentation> m_representation;
    std::shared_ptr<IFitness> m_fitness;
};
