//
// Created by konrad_guest on 08/01/2025.
// SMART


#include "benchmark/crossover_benchmark.h"
#include "utils/logging.h"

void CrossoverBenchmark::runBenchmark(const DNAInstance &instance)
{
    LOG_INFO("Starting crossover benchmark");
    
    // We create a list of crossovers to test
    std::vector<std::pair<std::string, std::shared_ptr<ICrossover>>> crossovers = {
        {"OnePoint",          std::make_shared<OnePointCrossover>()},
        {"OrderCrossover",    std::make_shared<OrderCrossover>()},
        {"EdgeRecombination", std::make_shared<EdgeRecombination>()},
        {"PMXCrossover",      std::make_shared<PMXCrossover>()}
    };

    // std::vector<std::pair<std::string, std::shared_ptr<ICrossover>>> crossovers = {
    //     {"OrderCrossover",    std::make_shared<OrderCrossover>()}
    // };

    std::cout << "\n=== Crossover Benchmark ===\n";
    std::cout << "Instance DNA length: " << instance.getDNA().size() << "\n";

    // For each crossover, run a GA once (or multiple times) 
    // and record the Levenshtein distance for the final solution.
    for (auto &co : crossovers) {
        const std::string &coName = co.first;
        auto coPtr = co.second;

        LOG_DEBUG("Running test case: " + coName);
        int distance = runOneGA(instance, coPtr, coName);

        std::cout << "Crossover: " << coName 
                  << ", Final Levenshtein distance = " << distance << "\n";
        LOG_INFO("Test results: " + std::to_string(distance));
    }
}

int CrossoverBenchmark::runOneGA(const DNAInstance& instance,
                                  std::shared_ptr<ICrossover> crossover,
                                  const std::string& outputFile)
{
    // Create a new config instance
    GAConfig config;
    if (!config.loadFromFile("config.cfg")) {
        std::cerr << "Failed to load GA configuration\n";
        return -1;
    }

    // Set up cache
    auto cache = std::make_shared<CachedPopulation>();
    config.setCache(cache);

    // Create GA components with config
    auto selection = std::make_shared<TournamentSelection>(config, cache);
    auto stopping = std::make_shared<MaxGenerationsStopping>(config);

    // Create and run GA
    GeneticAlgorithm ga(
        config.getRepresentation(),
        selection,
        crossover,
        config.getMutation(),
        config.getReplacement(),
        config.getFitness(),
        stopping,
        cache,
        config
    );

    ga.run(instance);

    // Save results if output file specified
    if (!outputFile.empty()) {
        std::ofstream out(outputFile);
        if (out.is_open()) {
            out << ga.getBestDNA() << std::endl;
            out.close();
        }
    }

    return static_cast<int>(ga.getBestFitness());
}
