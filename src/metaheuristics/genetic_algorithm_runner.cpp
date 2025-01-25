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

namespace {
    std::mutex configMutex;  // Declare mutex in anonymous namespace
    
    // Helper function to create a string from numeric value with proper error handling
    template<typename T>
    std::string safeToString(const T& value) {
        try {
            return std::to_string(value);
        } catch (const std::exception&) {
            return "error";
        }
    }
}

void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config)
{
    LOG_INFO("Updating configuration with instance parameters");
    
    try {
        // Store current algorithm parameters that we want to preserve from config.cfg
        LOG_DEBUG("Preserving algorithm parameters from config file");
        const int maxGen = config.getMaxGenerations();
        const int popSize = config.getPopulationSize();
        const double mutRate = config.getMutationRate();
        const double crossProb = config.getCrossoverProbability();
        const double replRatio = config.getReplacementRatio();
        
        LOG_DEBUG("Current algorithm parameters:");
        LOG_DEBUG("  Max generations: " + safeToString(maxGen));
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
        config.setMaxGenerations(maxGen);
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
                  ", crossoverType=" + config.getCrossoverType() +
                  ", selectionMethod=" + config.getSelectionMethod());
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
        // Load and validate configuration
        GAConfig config;
        {
            std::lock_guard<std::mutex> lock(configMutex);
            if (!config.loadFromFile(configFile)) {
                throw std::runtime_error("Failed to load configuration from file: " + configFile);
            }
        }
        
        // Update configuration with instance parameters
        updateConfigWithInstanceParams(instance, config);
        
        // Create and validate components with proper error handling
        std::shared_ptr<IPopulationCache> cache;
        std::shared_ptr<IRepresentation> representation;
        std::shared_ptr<ISelection> selection;
        std::shared_ptr<ICrossover> crossover;
        std::shared_ptr<IMutation> mutation;
        std::shared_ptr<IReplacement> replacement;
        std::shared_ptr<IFitness> fitness;
        std::shared_ptr<IStopping> stopping;
        
        try {
            cache = config.getCache();
            representation = config.getRepresentation();
            selection = config.getSelection();
            crossover = config.getCrossover(config.getCrossoverType());
            mutation = config.getMutation();
            replacement = config.getReplacement();
            fitness = config.getFitness();
            stopping = config.getStopping();
            
            // Validate all components
            if (!cache) throw std::runtime_error("Failed to create cache");
            if (!representation) throw std::runtime_error("Failed to create representation");
            if (!selection) throw std::runtime_error("Failed to create selection");
            if (!crossover) throw std::runtime_error("Failed to create crossover");
            if (!mutation) throw std::runtime_error("Failed to create mutation");
            if (!replacement) throw std::runtime_error("Failed to create replacement");
            if (!fitness) throw std::runtime_error("Failed to create fitness");
            if (!stopping) throw std::runtime_error("Failed to create stopping criteria");
            
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create GA components: " + std::string(e.what()));
        }
        
        // Create and run GA with proper error handling
        std::unique_ptr<GeneticAlgorithm> ga;
        try {
            ga = std::make_unique<GeneticAlgorithm>(
                representation, selection, crossover, mutation,
                replacement, fitness, stopping, cache, config
            );
            ga->setProcessId(processId);
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create genetic algorithm: " + std::string(e.what()));
        }
        
        // Run the algorithm
        try {
            ga->run(instance);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error during GA execution: " + std::string(e.what()));
        }
        
        // Save results with proper error handling
        try {
            const std::filesystem::path outputPath(outputFile);
            const auto outputDir = outputPath.parent_path();
            if (!outputDir.empty() && !std::filesystem::exists(outputDir)) {
                std::filesystem::create_directories(outputDir);
            }
            
            std::ofstream out(outputFile);
            if (!out.is_open()) {
                throw std::runtime_error("Cannot open output file for writing: " + outputFile);
            }
            
            out << std::fixed << std::setprecision(6);
            const auto bestFitness = ga->getBestFitness();
            const auto bestDNA = ga->getBestDNA();
            
            if (bestDNA.empty()) {
                throw std::runtime_error("Best DNA solution is empty");
            }
            
            out << "Best Fitness: " << bestFitness << "\n";
            out << "Best DNA: " << bestDNA << "\n";
            out.close();
            
            LOG_INFO("Process " + procId + ": Results saved successfully");
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to save results: " + std::string(e.what()));
        }
        
    } catch (const std::exception& e) {
        LOG_ERROR("Process " + procId + ": " + e.what());
        throw;
    }
    
    LOG_INFO("Process " + procId + ": Genetic algorithm execution completed");
} 