#include "metaheuristics/genetic_algorithm_runner.h"
#include "metaheuristics/genetic_algorithm.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "metaheuristics/population_cache.h"
#include "utils/logging.h"
#include <fstream>
#include <iostream>

void updateConfigWithInstanceParams(const DNAInstance& instance, GAConfig& config)
{
    // Store current algorithm parameters
    int currentMaxGenerations = config.getMaxGenerations();
    int currentPopulationSize = config.populationSize;
    double currentMutationRate = config.mutationRate;
    double currentReplacementRatio = config.replacementRatio;
    double currentCrossoverProbability = config.crossoverProbability;
    
    // Only update instance-specific parameters
    // DNA Generation parameters
    if (config.k <= 0) config.k = instance.getK();
    if (config.deltaK < 0) config.deltaK = instance.getDeltaK();
    if (config.lNeg < 0) config.lNeg = instance.getLNeg();
    if (config.lPoz < 0) config.lPoz = instance.getLPoz();
    if (!config.repAllowed) config.repAllowed = instance.isRepAllowed();
    if (config.probablePositive < 0) config.probablePositive = instance.getProbablePositive();
    
    // Restore algorithm parameters to ensure they weren't accidentally modified
    config.setMaxGenerations(currentMaxGenerations);
    config.populationSize = currentPopulationSize;
    config.mutationRate = currentMutationRate;
    config.replacementRatio = currentReplacementRatio;
    config.crossoverProbability = currentCrossoverProbability;
    
    // Log the current state of parameters
    LOG_DEBUG("Updated instance parameters: k=" + std::to_string(config.k) + 
              ", deltaK=" + std::to_string(config.deltaK) +
              ", lNeg=" + std::to_string(config.lNeg) +
              ", lPoz=" + std::to_string(config.lPoz));
              
    LOG_DEBUG("Preserved algorithm parameters: maxGenerations=" + std::to_string(config.getMaxGenerations()) +
              ", populationSize=" + std::to_string(config.populationSize) +
              ", mutationRate=" + std::to_string(config.mutationRate));
}

void runGeneticAlgorithm(const DNAInstance& instance, const std::string& outputFile, int processId, const std::string& difficulty)
{
    // Create a new config instance
    GAConfig config;
    if (!config.loadFromFile("config.cfg")) {
        std::cerr << "Failed to load GA configuration\n";
        return;
    }
    
    // Update instance-specific parameters
    updateConfigWithInstanceParams(instance, config);  // Pass config as parameter
    
    // Create and run GA with this config
    auto cache = std::make_shared<CachedPopulation>();
    config.setCache(cache);
    
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        config.getCrossover(config.crossoverType),
        config.getMutation(),
        config.getReplacement(),
        config.getFitness(),
        config.getStopping(),
        cache,
        config  // Pass config to GA
    );
    
    ga.run(instance);
    
    // Save results
    std::ofstream out(outputFile);
    if (out.is_open()) {
        out << ga.getBestDNA() << std::endl;
        out.close();
    }
} 