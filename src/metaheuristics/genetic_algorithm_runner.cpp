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
    LOG_DEBUG("Updating instance-specific parameters");

    // Only update instance parameters if they haven't been set in the config file
    // (i.e., they still have their default values)
    if (config.getK() == 7) {  // Default value
        config.setK(instance.getK());
    }
    if (config.getDeltaK() == 0) {  // Default value
        config.setDeltaK(instance.getDeltaK());
    }
    if (config.getLNeg() == 10) {  // Default value
        config.setLNeg(instance.getLNeg());
    }
    if (config.getLPoz() == 10) {  // Default value
        config.setLPoz(instance.getLPoz());
    }
    if (config.isRepAllowed() == true) {  // Default value
        config.setRepAllowed(instance.isRepAllowed());
    }
    if (config.getProbablePositive() == 0) {  // Default value
        config.setProbablePositive(instance.getProbablePositive());
    }

    LOG_INFO("Configuration updated successfully");
    LOG_DEBUG("Final configuration parameters:");
    LOG_DEBUG("  k=" + std::to_string(config.getK()) + 
              ", deltaK=" + std::to_string(config.getDeltaK()) + 
              ", lNeg=" + std::to_string(config.getLNeg()) + 
              ", lPoz=" + std::to_string(config.getLPoz()));
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
            LOG_INFO("Loading configuration from file: " + configFile);
            if (!config.loadFromFile(configFile)) {
                throw std::runtime_error("Failed to load configuration from file: " + configFile);
            }
            LOG_DEBUG("Configuration loaded successfully");
            LOG_DEBUG("Loaded parameters:");
            LOG_DEBUG("  populationSize=" + safeToString(config.getPopulationSize()) +
                     ", mutationRate=" + safeToString(config.getMutationRate()) +
                     ", crossoverProbability=" + safeToString(config.getCrossoverProbability()) +
                     ", replacementRatio=" + safeToString(config.getReplacementRatio()));
        }
        
        // Update only instance-specific parameters
        updateConfigWithInstanceParams(instance, config);
        
        LOG_DEBUG("Creating representation...");
        auto representation = std::make_unique<PermutationRepresentation>();
        if (!representation) {
            throw std::runtime_error("Failed to create representation");
        }
        
        LOG_DEBUG("Initializing cache...");
        auto cache = std::make_shared<PopulationCache>();
        if (!cache) {
            throw std::runtime_error("Failed to create cache");
        }
        cache->reserve(config.getPopulationSize() * 2);  // Reserve space for population and offspring
        cache->enableDiversityTracking(true);  // Enable diversity tracking
        cache->setDiversityThreshold(config.getDiversityParams().sharingRadius);
        config.setCache(cache);
        
        LOG_DEBUG("Cache initialized with diversity tracking enabled");
        
        LOG_DEBUG("Creating genetic algorithm...");
        GeneticAlgorithm ga(std::move(representation), config, debugMode);
        ga.setProcessId(processId);
        
        LOG_DEBUG("Starting genetic algorithm run...");
        std::string result = ga.run(instance);
        
        if (result.empty()) {
            LOG_WARNING("Genetic algorithm returned empty result");
        } else {
            LOG_INFO("Genetic algorithm completed with result length: " + std::to_string(result.length()));
        }
        
        // Save result if output file is specified
        if (!outputFile.empty()) {
            LOG_DEBUG("Saving result to file: " + outputFile);
            std::ofstream out(outputFile);
            if (out.is_open()) {
                out << result;
                out.close();
                LOG_INFO("Result saved successfully");
            } else {
                LOG_ERROR("Failed to open output file: " + outputFile);
            }
        }
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in genetic algorithm runner: " + std::string(e.what()));
        throw;
    } catch (...) {
        LOG_ERROR("Unknown error in genetic algorithm runner");
        throw;
    }
    
    LOG_INFO("Process " + procId + ": Genetic algorithm execution completed");
} 