//
// Created by konrad_guest on 5/01/2025.
// SMART


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
#include "metaheuristics/adaptive_crossover.h"
#include "benchmark/adaptive_crossover_benchmark.h"

void runAdaptiveCrossoverBenchmark(const DNAInstance& instance) {
    AdaptiveCrossoverBenchmark benchmark("adaptive_crossover_results.csv");
    benchmark.runBenchmark(instance);
}

void runGeneticAlgorithm(const DNAInstance& instance) {
    // Get the configuration
    auto& config = GAConfig::getInstance();
    config.setMutationRate(0.7);
    config.setMaxGenerations(10000);
    config.setPopulationSize(200);  // Keep it consistent
    config.setReplacementRatio(0.7);  // Keep 30% of parents
    
    // Create and set the cache
    auto cache = std::make_shared<CachedPopulation>();
    config.setCache(cache);  // This ensures cache is used in selection and replacement
    
    // Create GA with configuration
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        std::make_shared<AdaptiveCrossover>(),  // Use adaptive crossover
        config.getMutation(),
        config.getReplacement(),
        std::make_shared<OptimizedGraphBasedFitness>(),
        std::make_shared<MaxGenerationsStopping>(1200),
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
    // Opcjonalnie: seed rand(), jeśli jest gdzieś wykorzystywany
    srand(static_cast<unsigned int>(time(nullptr)));

    // Ustawiamy parametry zgodnie z nową logiką.
    // Przyjmijmy przykładowe wartości.
    // Jeśli wykraczają poza zakres, w DNAInstanceBuilder zostaną skorygowane do domyślnych.
    int n = 400;         // Długość DNA (300–700)
    int k = 8;           // Długość oligo (7–10)
    int deltaK = 1;      // Zakres zmian długości oligo (0–2)
    int lNeg = 10;       // Błędy negatywne (może być 0 lub >= 10)
    int lPoz = 10;       // Błędy pozytywne (może być 0 lub >= 10)
    bool repAllowed = true;   // Czy dozwolone powtórzenia
    int probablePos = 0;      // 0 - losowe sekwencje, 1 - duplikowanie + modyfikacje

    // Budujemy instancję przy pomocy buildera.
    DNAInstanceBuilder builder;
    builder.setN(n)
           .setK(k)
           .setDeltaK(deltaK)
           .setLNeg(lNeg)
           .setLPoz(lPoz)
           .setRepAllowed(repAllowed)
           .setProbablePositive(probablePos)

           // Generujemy DNA i spektrum
           .buildDNA()
           .buildSpectrum();

    // Wprowadzamy błędy negatywne, jeśli lNeg > 0
    // NegativeErrorIntroducer przyjmuje liczbę int jako lNeg
    if(lNeg > 0) {
        NegativeErrorIntroducer negErr(lNeg);
        builder.applyError(&negErr);
    }

    // Wprowadzamy błędy pozytywne, jeśli lPoz > 0
    // PositiveErrorIntroducer przyjmuje liczbę int jako lPoz
    if(lPoz > 0) {
        PositiveErrorIntroducer posErr(lPoz);
        builder.applyError(&posErr);
    }

    // Pobieramy gotową instancję
    DNAInstance instance = builder.getInstance();

    // Zapisujemy do pliku
    bool saved = InstanceIO::saveInstance(instance, "instance.txt");
    if(!saved) {
        std::cerr << "Bład zapisu instancji!\n";
    }
    else {
        std::cout << "Instancja zapisana do pliku 'instance.txt'.\n";
    }
    
    // Odczytujemy z pliku, by sprawdzić poprawność
    DNAInstance loadedInst;
    bool loaded = InstanceIO::loadInstance("instance.txt", loadedInst);
    if(!loaded) {
        std::cerr << "Bład odczytu instancji!\n";
        return 1;
    }
    
    std::cout << "Wczytano instancje z pliku. Parametry:\n"
              << "  n=" << loadedInst.getN()
              << ", k=" << loadedInst.getK()
              << ", deltaK=" << loadedInst.getDeltaK()
              << ", lNeg=" << loadedInst.getLNeg()
              << ", lPoz=" << loadedInst.getLPoz()
              << "\nDlugosc DNA=" << loadedInst.getDNA().size()
              << ", rozmiar spektrum=" << loadedInst.getSpectrum().size() 
              << std::endl;

    // Przykładowe wywołanie algorytmu genetycznego
    //runGeneticAlgorithm(loadedInst);
    //Profiler::getInstance().saveReport("profiling_report.csv");

    runAdaptiveCrossoverBenchmark(loadedInst);

    //  NaiveBenchmark nb;
    //  nb.runBenchmark(loadedInst);
    //
    // CrossoverBenchmark cb;
    // cb.runBenchmark(loadedInst);

    return 0;
}
