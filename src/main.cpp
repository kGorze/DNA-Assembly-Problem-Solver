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



int main()
{
    DNAInstanceBuilder builder;
    builder.setN(500).setK(9)
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

    // 2) Possibly run naive benchmark
    std::cout << "\n=== Uruchamiam NaiveBenchmark ===\n";
    NaiveBenchmark nb;
    nb.runBenchmark(loadedInst);

    // 3) Setup metaheuristic with representation + strategies
    auto representation = std::make_shared<DirectDNARepresentation>();
    auto selection  = std::make_shared<TournamentSelection>(3);
    auto crossover  = std::make_shared<OnePointCrossover>();
    auto mutation   = std::make_shared<PointMutation>();
    auto replacement = std::make_shared<FullReplacement>();
    auto fitness    = std::make_shared<SimpleFitness>();
    auto stopping   = std::make_shared<MaxGenerationsStopping>(50);

    // 4) Create GeneticAlgorithm
    GeneticAlgorithm ga(
        representation, // 7 args total
        selection,
        crossover,
        mutation,
        replacement,
        fitness,
        stopping
    );

    // 5) Run
    ga.run(loadedInst);

    return 0;
}

