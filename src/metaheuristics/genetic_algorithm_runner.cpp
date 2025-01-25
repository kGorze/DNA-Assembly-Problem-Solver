#include "../include/metaheuristics/genetic_algorithm_runner.h"
#include "../include/metaheuristics/genetic_algorithm.h"
#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../include/metaheuristics/population_cache_impl.h"
#include "../include/utils/logging.h"
#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <mutex>

void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config)
{
    LOG_INFO("Updating configuration with instance parameters");
    
    // Store current algorithm parameters that we want to preserve from config.cfg
    LOG_DEBUG("Preserving algorithm parameters from config file");
    int maxGen = config.getMaxGenerations();
    int popSize = config.getPopulationSize();
    double mutRate = config.getMutationRate();
    double crossProb = config.getCrossoverProbability();
    double replRatio = config.getReplacementRatio();
    
    LOG_DEBUG("Current algorithm parameters:");
    LOG_DEBUG("  Max generations: " + std::to_string(maxGen));
    LOG_DEBUG("  Population size: " + std::to_string(popSize));
    LOG_DEBUG("  Mutation rate: " + std::to_string(mutRate));
    LOG_DEBUG("  Crossover probability: " + std::to_string(crossProb));
    LOG_DEBUG("  Replacement ratio: " + std::to_string(replRatio));
    
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
    config.setMaxGenerations(maxGen);
    config.setPopulationSize(popSize);
    config.setMutationRate(mutRate);
    config.setCrossoverProbability(crossProb);
    config.setReplacementRatio(replRatio);
    
    LOG_INFO("Configuration updated successfully");
    LOG_DEBUG("Final configuration parameters:");
    LOG_DEBUG("  k=" + std::to_string(config.getK()) +
              ", deltaK=" + std::to_string(config.getDeltaK()) +
              ", lNeg=" + std::to_string(config.getLNeg()) +
              ", lPoz=" + std::to_string(config.getLPoz()));
    LOG_DEBUG("  populationSize=" + std::to_string(config.getPopulationSize()) +
              ", mutationRate=" + std::to_string(config.getMutationRate()) +
              ", crossoverType=" + config.getCrossoverType() +
              ", selectionMethod=" + config.getSelectionMethod());
}

void runGeneticAlgorithm(
    const DNAInstance& instance,
    const std::string& outputFile,
    int processId,
    const std::string& configFile,
    bool debugMode)
{
    LOG_INFO("Process " + std::to_string(processId) + ": Starting genetic algorithm execution");
    LOG_INFO("Process " + std::to_string(processId) + ": Configuration file: " + configFile);
    LOG_INFO("Process " + std::to_string(processId) + ": Output file: " + outputFile);
    LOG_INFO("Process " + std::to_string(processId) + ": Debug mode: " + std::string(debugMode ? "enabled" : "disabled"));
    
    try {
        // Validate paths
        if (configFile.empty()) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Empty config file path");
            throw std::runtime_error("Empty config file path");
        }
        
        // Check if config file exists
        std::ifstream configCheck(configFile);
        if (!configCheck.good()) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Config file does not exist or is not accessible: " + configFile);
            throw std::runtime_error("Config file not accessible: " + configFile);
        }
        configCheck.close();
        
        // Check if output directory exists and is writable
        std::filesystem::path outputPath(outputFile);
        auto outputDir = outputPath.parent_path();
        if (!outputDir.empty() && !std::filesystem::exists(outputDir)) {
            LOG_INFO("Process " + std::to_string(processId) + ": Creating output directory: " + outputDir.string());
            std::filesystem::create_directories(outputDir);
        }
        
        // Try to create a test file in output directory
        std::ofstream testFile(outputFile + ".test");
        if (!testFile.is_open()) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Cannot write to output directory: " + outputDir.string());
            throw std::runtime_error("Cannot write to output directory");
        }
        testFile.close();
        std::filesystem::remove(outputFile + ".test");
        
        // Load configuration
        GAConfig config;
        LOG_INFO("Process " + std::to_string(processId) + ": Loading configuration from " + configFile);
        try {
            if (!config.loadFromFile(configFile)) {
                LOG_ERROR("Process " + std::to_string(processId) + ": Failed to load configuration from " + configFile);
                throw std::runtime_error("Configuration loading failed");
            }
            LOG_INFO("Process " + std::to_string(processId) + ": Configuration loaded successfully");
            LOG_INFO("Process " + std::to_string(processId) + ": Configuration parameters:");
            LOG_INFO("- Population size: " + std::to_string(config.getPopulationSize()));
            LOG_INFO("- Max generations: " + std::to_string(config.getMaxGenerations()));
            LOG_INFO("- Mutation rate: " + std::to_string(config.getMutationRate()));
            LOG_INFO("- Crossover probability: " + std::to_string(config.getCrossoverProbability()));
            LOG_INFO("- Selection method: " + config.getSelectionMethod());
            LOG_INFO("- Crossover type: " + config.getCrossoverType());
        } catch (const std::exception& e) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Exception while loading config: " + std::string(e.what()));
            throw;
        }

        // Validate instance
        if (instance.getSpectrum().empty()) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Empty spectrum in instance");
            throw std::runtime_error("Invalid instance - empty spectrum");
        }
        
        // Update config with instance parameters
        LOG_INFO("Process " + std::to_string(processId) + ": Updating configuration with instance parameters");
        try {
            updateConfigWithInstanceParams(instance, config);
            LOG_INFO("Process " + std::to_string(processId) + ": Updated instance parameters:");
            LOG_INFO("- K: " + std::to_string(config.getK()));
            LOG_INFO("- Delta K: " + std::to_string(config.getDeltaK()));
            LOG_INFO("- L Neg: " + std::to_string(config.getLNeg()));
            LOG_INFO("- L Poz: " + std::to_string(config.getLPoz()));
        } catch (const std::exception& e) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Failed to update config with instance params: " + std::string(e.what()));
            throw;
        }

        // Create components
        LOG_INFO("Process " + std::to_string(processId) + ": Creating genetic algorithm components");
        
        auto cache = config.getCache();
        LOG_INFO("Process " + std::to_string(processId) + ": Created population cache");
        
        auto representation = config.getRepresentation();
        LOG_INFO("Process " + std::to_string(processId) + ": Created representation");
        
        auto selection = config.getSelection();
        LOG_INFO("Process " + std::to_string(processId) + ": Created selection operator");
        
        auto crossover = config.getCrossover(config.getCrossoverType());
        LOG_INFO("Process " + std::to_string(processId) + ": Created crossover operator");
        
        auto mutation = config.getMutation();
        LOG_INFO("Process " + std::to_string(processId) + ": Created mutation operator");
        
        auto replacement = config.getReplacement();
        LOG_INFO("Process " + std::to_string(processId) + ": Created replacement operator");
        
        auto fitness = config.getFitness();
        LOG_INFO("Process " + std::to_string(processId) + ": Created fitness evaluator");
        
        auto stopping = config.getStopping();
        LOG_INFO("Process " + std::to_string(processId) + ": Created stopping criteria");

        // Create and run GA
        LOG_INFO("Process " + std::to_string(processId) + ": Creating genetic algorithm instance");
        GeneticAlgorithm ga(
            representation,
            selection,
            crossover,
            mutation,
            replacement,
            fitness,
            stopping,
            cache,
            config
        );
        LOG_INFO("Process " + std::to_string(processId) + ": Genetic algorithm instance created");

        LOG_INFO("Process " + std::to_string(processId) + ": Starting genetic algorithm run");
        ga.run(instance);
        LOG_INFO("Process " + std::to_string(processId) + ": Genetic algorithm run completed");

        // Save results
        LOG_INFO("Process " + std::to_string(processId) + ": Saving results to " + outputFile);
        std::ofstream out(outputFile);
        if (!out.is_open()) {
            LOG_ERROR("Process " + std::to_string(processId) + ": Failed to open output file for writing");
            throw std::runtime_error("Cannot open output file for writing");
        }
        
        out << "Best Fitness: " << ga.getBestFitness() << "\n";
        out << "Best DNA: " << ga.getBestDNA() << "\n";
        out.close();
        LOG_INFO("Process " + std::to_string(processId) + ": Results saved successfully");

    } catch (const std::exception& e) {
        LOG_ERROR("Process " + std::to_string(processId) + ": Exception in runGeneticAlgorithm: " + std::string(e.what()));
        throw;
    }
    LOG_INFO("Process " + std::to_string(processId) + ": Genetic algorithm execution completed");
} 