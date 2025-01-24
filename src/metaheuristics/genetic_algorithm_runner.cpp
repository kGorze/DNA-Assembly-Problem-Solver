#include "../include/metaheuristics/genetic_algorithm_runner.h"
#include "../include/metaheuristics/genetic_algorithm.h"
#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../include/metaheuristics/population_cache_impl.h"
#include "../include/utils/logging.h"
#include <fstream>
#include <iostream>
#include <thread>
#include <chrono>

void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config)
{
    // Store current algorithm parameters
    int currentPopulationSize = config.getPopulationSize();
    double currentMutationRate = config.getMutationRate();
    double currentCrossoverProbability = config.getCrossoverProbability();
    
    // Only update instance-specific parameters
    if (config.getK() <= 0) config.setK(instance.getK());
    if (config.getDeltaK() < 0) config.setDeltaK(instance.getDeltaK());
    if (config.getLNeg() < 0) config.setLNeg(instance.getLNeg());
    if (config.getLPoz() < 0) config.setLPoz(instance.getLPoz());
    if (!config.getRepAllowed()) config.setRepAllowed(instance.isRepAllowed());
    if (config.getProbablePositive() < 0) config.setProbablePositive(instance.getProbablePositive());
    
    // Restore algorithm parameters
    config.setPopulationSize(currentPopulationSize);
    config.setMutationRate(currentMutationRate);
    config.setCrossoverProbability(currentCrossoverProbability);
    
    std::string debugMsg = "Updated instance parameters: k=" + std::to_string(config.getK()) +
                          ", deltaK=" + std::to_string(config.getDeltaK()) +
                          ", lNeg=" + std::to_string(config.getLNeg()) +
                          ", lPoz=" + std::to_string(config.getLPoz());
    LOG_DEBUG(debugMsg);
              
    std::string restoreMsg = "Restored algorithm parameters: populationSize=" + 
                            std::to_string(config.getPopulationSize()) +
                            ", mutationRate=" + std::to_string(config.getMutationRate());
    LOG_DEBUG(restoreMsg);
}

void runGeneticAlgorithm(const DNAInstance& instance, const std::string& outputFile, int seed, const std::string& configFile)
{
    // Get singleton instance
    auto& config = GAConfig::getInstance();
    
    try {
        config.loadFromFile(configFile);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load GA configuration: " << e.what() << std::endl;
        return;
    }
    
    // Update instance-specific parameters
    updateConfigWithInstanceParams(instance, config);
    
    // Create and run GA with this config
    auto cache = std::make_shared<SimplePopulationCache>();
    config.setCache(cache);
    
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        config.getCrossover(config.getCrossoverType()),
        config.getMutation(),
        config.getReplacement(),
        config.getFitness(),
        config.getStopping(),
        cache,
        config
    );
    
    ga.run(instance);
    
    // Save results
    std::ofstream out(outputFile);
    if (out.is_open()) {
        out << ga.getBestDNA() << std::endl;
        out.close();
    }
} 