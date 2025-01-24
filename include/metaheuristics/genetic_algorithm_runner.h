#ifndef GENETIC_ALGORITHM_RUNNER_H
#define GENETIC_ALGORITHM_RUNNER_H

#include <string>
#include "dna/dna_instance.h"
#include "configuration/genetic_algorithm_configuration.h"

// Function to update config with instance-specific parameters
void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config);

// Function to run the genetic algorithm with given parameters
void runGeneticAlgorithm(const DNAInstance& instance, 
                        const std::string& outputFile, 
                        int processId, 
                        const std::string& difficulty);

#endif // GENETIC_ALGORITHM_RUNNER_H 