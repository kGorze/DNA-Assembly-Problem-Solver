//
// Created by konrad_guest on 08/01/2025.
//

#ifndef CROSSOVER_BENCHMARK_H
#define CROSSOVER_BENCHMARK_H

#include "metaheuristics/representation.h"
#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"
#include "metaheuristics/genetic_algorithm.h"
#include "naive/naive_reconstruction.h"
#include "configuration/genetic_algorithm_configuration.h"


#include <string>
#include <vector>
#include <memory>
#include <iostream>

class CrossoverBenchmark {
public:
    void runBenchmark(const DNAInstance &instance);

private:
    // Helper function to run one GA with a given crossover
    int runOneGA(const DNAInstance &instance, 
                 std::shared_ptr<ICrossover> crossover, 
                 const std::string &crossoverName);
};


#endif //CROSSOVER_BENCHMARK_H
