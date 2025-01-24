//
// Created by konrad_guest on 08/01/2025.
// SMART

#ifndef CROSSOVER_BENCHMARK_H
#define CROSSOVER_BENCHMARK_H

#include "../interfaces/i_representation.h"
#include "../interfaces/i_selection.h"
#include "../interfaces/i_crossover.h"
#include "../interfaces/i_mutation.h"
#include "../interfaces/i_replacement.h"
#include "../interfaces/i_fitness.h"
#include "../interfaces/i_stopping.h"
#include "../interfaces/i_population_cache.h"

#include "../metaheuristics/representation_impl.h"
#include "../metaheuristics/selection_impl.h"
#include "../metaheuristics/crossover_impl.h"
#include "../metaheuristics/mutation_impl.h"
#include "../metaheuristics/replacement_impl.h"
#include "../metaheuristics/fitness_impl.h"
#include "../metaheuristics/stopping_criteria_impl.h"
#include "../metaheuristics/genetic_algorithm.h"
#include "../naive/naive_reconstruction.h"
#include "../configuration/genetic_algorithm_configuration.h"

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
