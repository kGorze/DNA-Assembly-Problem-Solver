#include "../include/metaheuristics/genetic_algorithm_runner.h"
#include "../include/metaheuristics/genetic_algorithm.h"
#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../include/metaheuristics/population_cache_impl.h"
#include "../include/utils/logging.h"
#include "../include/metaheuristics/representation.h"
#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <mutex>

namespace {
    // Helper function to safely convert numbers to strings
    template<typename T>
    std::string safeToString(T value) {
        try {
            return std::to_string(value);
        } catch (const std::exception& e) {
            LOG_ERROR("Error converting value to string: " + std::string(e.what()));
            return "";
        }
    }
}

void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config)
{
    LOG_INFO("Updating configuration with instance parameters");
    
    try {
        // Store current algorithm parameters that we want to preserve from config.cfg
        LOG_DEBUG("Preserving algorithm parameters from config file");
        const int popSize = config.getPopulationSize();
        const double mutRate = config.getMutationRate();
        const double crossProb = config.getCrossoverProbability();
        const double replRatio = config.getReplacementRatio();
        
        LOG_DEBUG("Current algorithm parameters:");
        LOG_DEBUG("  Population size: " + safeToString(popSize));
        LOG_DEBUG("  Mutation rate: " + safeToString(mutRate));
        LOG_DEBUG("  Crossover probability: " + safeToString(crossProb));
        LOG_DEBUG("  Replacement ratio: " + safeToString(replRatio));
        
        // Update all instance-specific parameters from the instance
        LOG_DEBUG("Updating instance-specific parameters");
        config.setK(instance.getK());
        config.setDeltaK(instance.getDeltaK());
        config.setLNeg(instance.getLNeg());
        config.setLPoz(instance.getLPoz());
        config.setRepAllowed(instance.isRepAllowed());
        config.setProbablePositive(instance.getProbablePositive());
        
        // Restore algorithm parameters from config.cfg
        LOG_DEBUG("Restoring algorithm parameters");
        config.setPopulationSize(popSize);
        config.setMutationRate(mutRate);
        config.setCrossoverProbability(crossProb);
        config.setReplacementRatio(replRatio);
        
        LOG_INFO("Configuration updated successfully");
        LOG_DEBUG("Final configuration parameters:");
        LOG_DEBUG("  k=" + safeToString(config.getK()) +
                  ", deltaK=" + safeToString(config.getDeltaK()) +
                  ", lNeg=" + safeToString(config.getLNeg()) +
                  ", lPoz=" + safeToString(config.getLPoz()));
        LOG_DEBUG("  populationSize=" + safeToString(config.getPopulationSize()) +
                  ", mutationRate=" + safeToString(config.getMutationRate()) +
                  ", crossoverProbability=" + safeToString(config.getCrossoverProbability()) +
                  ", replacementRatio=" + safeToString(config.getReplacementRatio()));
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to update configuration parameters: " + std::string(e.what()));
        throw;
    }
}

void runGeneticAlgorithm(
    const DNAInstance& instance,
    const std::string& outputFile,
    int processId,
    const std::string& configFile,
    bool debugMode)
{
    const std::string procId = safeToString(processId);
    LOG_INFO("Process " + procId + ": Starting genetic algorithm execution");
    
    try {
        // Load configuration
        GAConfig config;
        if (!configFile.empty()) {
            if (!config.loadFromFile(configFile)) {
                throw std::runtime_error("Failed to load configuration from file: " + configFile);
            }
        }
        
        // Update configuration with instance parameters
        updateConfigWithInstanceParams(instance, config);
        
        // Create representation
        auto representation = std::make_unique<PermutationRepresentation>();
        
        // Create and run genetic algorithm
        GeneticAlgorithm ga(std::move(representation), config, debugMode);
        ga.setProcessId(processId);
        
        // Run the algorithm
        std::string result = ga.run(instance);
        
        // Save result if output file is specified
        if (!outputFile.empty()) {
            std::ofstream out(outputFile);
            if (out.is_open()) {
                out << result;
                out.close();
            } else {
                LOG_ERROR("Failed to open output file: " + outputFile);
            }
        }
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in genetic algorithm runner: " + std::string(e.what()));
        throw;
    }
    
    LOG_INFO("Process " + procId + ": Genetic algorithm execution completed");
} 