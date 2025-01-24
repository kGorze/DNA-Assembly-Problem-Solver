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

void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config)
{
    // Store current algorithm parameters that we want to preserve from config.cfg
    int maxGen = config.getMaxGenerations();
    int popSize = config.getPopulationSize();
    double mutRate = config.getMutationRate();
    double crossProb = config.getCrossoverProbability();
    double replRatio = config.getReplacementRatio();
    
    // Update all instance-specific parameters from the instance
    config.setK(instance.getK());
    config.setDeltaK(instance.getDeltaK());
    config.setLNeg(instance.getLNeg());
    config.setLPoz(instance.getLPoz());
    config.setRepAllowed(instance.isRepAllowed());
    config.setProbablePositive(instance.getProbablePositive());
    
    // Restore algorithm parameters from config.cfg
    config.setMaxGenerations(maxGen);
    config.setPopulationSize(popSize);
    config.setMutationRate(mutRate);
    config.setCrossoverProbability(crossProb);
    config.setReplacementRatio(replRatio);
    
    LOG_DEBUG("Updated instance parameters from instance file:");
    LOG_DEBUG("  k=" + std::to_string(config.getK()) +
              ", deltaK=" + std::to_string(config.getDeltaK()) +
              ", lNeg=" + std::to_string(config.getLNeg()) +
              ", lPoz=" + std::to_string(config.getLPoz()));
              
    LOG_DEBUG("Preserved algorithm parameters from config.cfg:");
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
    if (debugMode) {
        Logger::setLogLevel(LogLevel::DEBUG);
        LOG_INFO("Debug mode enabled - detailed logging will be shown");
        LOG_INFO("==============================");
        LOG_INFO("Configuration:");
        LOG_INFO("  Config file: " + configFile);
        LOG_INFO("  Output file: " + outputFile);
        LOG_INFO("  Process ID: " + std::to_string(processId));
    }

    try {
        // Load configuration
        GAConfig& config = GAConfig::getInstance();
        if (!config.loadFromFile(configFile)) {
            LOG_ERROR("Failed to load configuration from " + configFile);
            throw std::runtime_error("Configuration loading failed");
        }

        if (debugMode) {
            LOG_INFO("Configuration loaded successfully");
            LOG_INFO("Algorithm parameters:");
            LOG_INFO("  Population size: " + std::to_string(config.getPopulationSize()));
            LOG_INFO("  Max generations: " + std::to_string(config.getMaxGenerations()));
            LOG_INFO("  Mutation rate: " + std::to_string(config.getMutationRate()));
            LOG_INFO("  Crossover probability: " + std::to_string(config.getCrossoverProbability()));
            LOG_INFO("  Replacement ratio: " + std::to_string(config.getReplacementRatio()));
            LOG_INFO("  Selection method: " + config.getSelectionType());
            LOG_INFO("  Crossover type: " + config.getCrossoverType());
            LOG_INFO("  Mutation method: " + config.getMutationType());
            LOG_INFO("");
            LOG_INFO("Instance parameters:");
            LOG_INFO("  N: " + std::to_string(instance.getN()));
            LOG_INFO("  K: " + std::to_string(instance.getK()));
            LOG_INFO("  Delta K: " + std::to_string(instance.getDeltaK()));
            LOG_INFO("  L Neg: " + std::to_string(instance.getLNeg()));
            LOG_INFO("  L Poz: " + std::to_string(instance.getLPoz()));
            LOG_INFO("  Spectrum size: " + std::to_string(instance.getSpectrum().size()));
            LOG_INFO("  Start index: " + std::to_string(instance.getStartIndex()));
            LOG_INFO("  Original DNA: " + instance.getDNA());
            LOG_INFO("==============================");
        }

        // Create population cache
        auto cache = std::make_shared<SimplePopulationCache>();
        config.setCache(cache);
        if (debugMode) LOG_INFO("Population cache initialized");

        // Create and initialize the genetic algorithm with all required components
        if (debugMode) LOG_INFO("Initializing genetic algorithm components...");
        auto representation = config.getRepresentation();
        auto selection = config.getSelection();
        auto crossover = config.getCrossover(config.getCrossoverType());
        auto mutation = config.getMutation();
        auto replacement = config.getReplacement();
        auto fitness = config.getFitness();
        auto stopping = config.getStopping();
        
        if (debugMode) {
            LOG_INFO("Components created:");
            LOG_INFO("  Representation: " + std::string(typeid(*representation).name()));
            LOG_INFO("  Selection: " + std::string(typeid(*selection).name()));
            LOG_INFO("  Crossover: " + std::string(typeid(*crossover).name()));
            LOG_INFO("  Mutation: " + std::string(typeid(*mutation).name()));
            LOG_INFO("  Replacement: " + std::string(typeid(*replacement).name()));
            LOG_INFO("  Fitness: " + std::string(typeid(*fitness).name()));
            LOG_INFO("  Stopping: " + std::string(typeid(*stopping).name()));
        }
        
        if (debugMode) LOG_INFO("Creating genetic algorithm instance...");
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
        ga.setProcessId(processId);
        if (debugMode) LOG_INFO("Genetic algorithm initialized successfully");
        
        // Open output file before starting
        std::ofstream out(outputFile);
        if (!out.is_open()) {
            throw std::runtime_error("Failed to open output file: " + outputFile);
        }
        if (debugMode) LOG_INFO("Output file opened: " + outputFile);

        // Start timing
        auto start = std::chrono::high_resolution_clock::now();
        if (debugMode) LOG_INFO("Starting genetic algorithm execution...");

        // Run the algorithm with the instance
        ga.run(instance);
        
        // End timing
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        if (debugMode) {
            // Log final results
            LOG_INFO("=== Final Results ===");
            LOG_INFO("Execution time: " + std::to_string(duration.count()) + "ms");
            LOG_INFO("Best fitness: " + std::to_string(ga.getBestFitness()));
            LOG_INFO("Best DNA: " + ga.getBestDNA());
            
            // Save results to output file
            out << "Execution time (ms): " << duration.count() << std::endl;
            out << "Best fitness: " << ga.getBestFitness() << std::endl;
            out << "Best DNA: " << ga.getBestDNA() << std::endl;
            
            LOG_INFO("Results saved to " + outputFile);
        }
        
        out.close();
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in genetic algorithm execution: " + std::string(e.what()));
        throw;
    }
} 