#pragma once

#include <string>
#include "dna/dna_instance.h"
#include "configuration/genetic_algorithm_configuration.h"

/*
    Updates the configuration with instance-specific parameters while preserving algorithm parameters
*/
void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config);

/*
    Runs the genetic algorithm with the given instance and configuration
    @param instance The DNA instance to solve
    @param outputFile Path to save results
    @param processId Process identifier (for parallel execution)
    @param configFile Path to configuration file
    @param debugMode Whether to enable detailed debug output
*/
void runGeneticAlgorithm(const DNAInstance& instance,
                        const std::string& outputFile,
                        int processId,
                        const std::string& configFile,
                        bool debugMode = false); 