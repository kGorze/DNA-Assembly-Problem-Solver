#include <iostream>
#include <random>
#include <chrono>

#include "generator/dna_generator.h"
#include "naive/naive_reconstruction.h"
#include "benchmark/naive_benchmark.h"

#include "metaheuristics/crossover.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/genetic_algorithm.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/selection.h"
#include "metaheuristics/stopping_criteria.h"

#include "benchmark/crossover_benchmark.h"

#include "configuration/genetic_algorithm_configuration.h"


void runGeneticAlgorithm(const DNAInstance& instance) {
    // Get the configuration
    auto& config = GAConfig::getInstance();
    config.setMutationRate(0.15);
    config.setPopulationSize(100);
    config.setMaxGenerations(197);
    config.setPopulationSize(100);  // Keep it consistent
    config.setReplacementRatio(0.7);  // Keep 30% of parents
    

    
    // Create GA with configuration
    auto cache = std::make_shared<CachedPopulation>();
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        config.getCrossover("order"),
        config.getMutation(),
        config.getReplacement(),
        std::make_shared<SmithWatermanFitness>(),
        std::make_shared<NoImprovementStopping>(300),
        cache
    );
    
    // Run the GA
    ga.run(instance);
    
    // Get and print results
    std::string reconstructedDNA = ga.getBestDNA();
    std::string originalDNA = instance.getDNA();
    
    std::cout << "\nOriginal DNA (first 100 bases): " 
              << originalDNA.substr(0, 100) << "...\n";
    std::cout << "Reconstructed DNA (first 100 bases): " 
              << reconstructedDNA.substr(0, 100) << "...\n";
    
    int distance = levenshteinDistance(originalDNA, reconstructedDNA);
    std::cout << "Levenshtein distance: " << distance << "\n";
}


int main()
{
    srand(static_cast<unsigned int>(time(nullptr))); // Seed rand()
    
    DNAInstanceBuilder builder;
    builder.setN(50).setK(4)
           .buildDNA()
           .buildSpectrum();
    
    NegativeErrorIntroducer negErr(5.0);
    builder.applyError(&negErr);
    
    PositiveErrorIntroducer posErr(3.0, 0.7);
    builder.applyError(&posErr);
    
    DNAInstance instance = builder.getInstance();
    bool saved = InstanceIO::saveInstance(instance, "instance.txt");
    if(!saved) {
        std::cerr << "Błąd zapisu!\n";
    }
    
    DNAInstance loadedInst;
    bool loaded = InstanceIO::loadInstance("instance.txt", loadedInst);
    if(loaded) {
        std::cout << "Wczytano instancje. n=" << loadedInst.getN()
                  << ", k=" << loadedInst.getK()
                  << ", dlugosc DNA=" << loadedInst.getDNA().size()
                  << ", rozmiar spektrum=" << loadedInst.getSpectrum().size() 
                  << std::endl;
    }
    
    // NaiveBenchmark nb;
    // nb.runBenchmark(loadedInst);
    
    //CrossoverBenchmark cb;
    //cb.runBenchmark(loadedInst);

    runGeneticAlgorithm(loadedInst);
    Profiler::getInstance().saveReport("profiling_report.csv");

    



    return 0;
}

