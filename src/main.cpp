#include <iostream> 
#include <random>
#include <chrono>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>

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

void runGeneticAlgorithm(const DNAInstance& instance, const std::string& outputFile = "", int processId = 0);

void printUsage() {
    std::cout << "Usage: dna_reconstruction <mode> [options]\n\n"
              << "Modes:\n"
              << "  debug             - Run in debug mode with default settings\n"
              << "  generate_instance - Generate test instances\n"
              << "  test_instance     - Solve DNA reconstruction from input file\n\n"
              << "Options for generate_instance:\n"
              << "  -n <value>        - DNA length (300-700)\n"
              << "  -k <value>        - Oligo length (7-10)\n"
              << "  -dk <value>       - Delta K range (0-2)\n"
              << "  -ln <value>       - Negative errors (0 or >=10)\n"
              << "  -lp <value>       - Positive errors (0 or >=10)\n"
              << "  -o <filename>     - Output filename\n\n"
              << "Options for test_instance:\n"
              << "  -i <filename>     - Input instance file\n"
              << "  -o <filename>     - Output results file\n"
              << "  -pid <value>      - Process ID (unique identifier)\n";
}

void runGeneticAlgorithm(const DNAInstance& instance, const std::string& outputFile, int processId) {
    std::cerr << "Starting GA with processId: " << processId << std::endl;
    auto& config = GAConfig::getInstance();
    config.setMutationRate(0.7);
    config.setMaxGenerations(100);
    config.setPopulationSize(100);
    config.setReplacementRatio(0.7);
    
    auto cache = std::make_shared<CachedPopulation>();
    config.setCache(cache);
    
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        std::make_shared<AdaptiveCrossover>(),
        config.getMutation(),
        config.getReplacement(),
        std::make_shared<OptimizedGraphBasedFitness>(),
        config.getStopping(),
        cache
    );
    
    
    // Ustawienie callbacku do aktualizacji postępu
    ga.setProgressCallback([processId](int generation, int maxGenerations, double bestFitness) {
        // Format status message do przekazania
        std::ostringstream status;
        status << std::fixed << std::setprecision(2) 
               << "Gen " << generation << "/" << maxGenerations 
               << " Fit: " << bestFitness;
               
        // Wypisanie komunikatu postępu
        std::cout << "PROGRESS_UPDATE:" << processId << ":" 
                  << ((double)generation / maxGenerations * 100.0) << ":" 
                  << status.str() << ":" << bestFitness << std::endl;
        std::cout.flush();
    });

    // Ustawienie identyfikatora procesu
    ga.setProcessId(processId); 
        
    ga.run(instance);
    
    std::string reconstructedDNA = ga.getBestDNA();
    std::string originalDNA = instance.getDNA();
    
    int distance = levenshteinDistance(originalDNA, reconstructedDNA);
    
    if (outputFile.empty()) {
        std::cout << "\nOriginal DNA (first 100 bases): " 
                  << originalDNA.substr(0, 100) << "...\n";
        std::cout << "Reconstructed DNA (first 100 bases): " 
                  << reconstructedDNA.substr(0, 100) << "...\n";
        std::cout << "Levenshtein distance: " << distance << "\n";
    } else {
        std::ofstream outFile(outputFile);
        if (outFile.is_open()) {
            outFile << "Original DNA: " << originalDNA << "\n";
            outFile << "Reconstructed DNA: " << reconstructedDNA << "\n";
            outFile << "Levenshtein distance: " << distance << "\n";
            outFile.close();
            
            // Wypisanie finalnego statusu
            std::cout << "PROGRESS_UPDATE:" << processId << ":100:Completed:0\n";
            std::cout.flush();
        } else {
            std::cerr << "Failed to open output file: " << outputFile << std::endl;
        }
    }
}

bool generateInstance(int n, int k, int deltaK, int lNeg, int lPoz, const std::string& outputFile) {
    DNAInstanceBuilder builder;
    builder.setN(n)
           .setK(k)
           .setDeltaK(deltaK)
           .setLNeg(lNeg)
           .setLPoz(lPoz)
           .setRepAllowed(true)
           .setProbablePositive(0)
           .buildDNA()
           .buildSpectrum();

    if(lNeg > 0) {
        NegativeErrorIntroducer negErr(lNeg);
        builder.applyError(&negErr);
    }

    if(lPoz > 0) {
        PositiveErrorIntroducer posErr(lPoz);
        builder.applyError(&posErr);
    }

    DNAInstance instance = builder.getInstance();
    return InstanceIO::saveInstance(instance, outputFile);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printUsage();
        return 1;
    }

    std::string mode = argv[1];

    if (mode == "debug") {
        // Domyślne parametry dla trybu debug
        int n = 400;
        int k = 8;
        int deltaK = 1;
        int lNeg = 10;
        int lPoz = 10;
        
        DNAInstanceBuilder builder;
        builder.setN(n)
               .setK(k)
               .setDeltaK(deltaK)
               .setLNeg(lNeg)
               .setLPoz(lPoz)
               .setRepAllowed(true)
               .setProbablePositive(0)
               .buildDNA()
               .buildSpectrum();

        if(lNeg > 0) {
            NegativeErrorIntroducer negErr(lNeg);
            builder.applyError(&negErr);
        }

        if(lPoz > 0) {
            PositiveErrorIntroducer posErr(lPoz);
            builder.applyError(&posErr);
        }

        DNAInstance instance = builder.getInstance();
        bool saved = InstanceIO::saveInstance(instance, "debug_instance.txt");
        
        if(!saved) {
            std::cerr << "Failed to save debug instance!\n";
            return 1;
        }

        runGeneticAlgorithm(instance, "");
        Profiler::getInstance().saveReport("debug_profiling_report.csv");
        
    } else if (mode == "generate_instance") {
        int n = 400, k = 8, deltaK = 1, lNeg = 10, lPoz = 10;
        std::string outputFile = "generated_instance.txt";
        
        for (int i = 2; i < argc; i += 2) {
            if (i + 1 >= argc) break;
            
            if (strcmp(argv[i], "-n") == 0) n = std::stoi(argv[i+1]);
            else if (strcmp(argv[i], "-k") == 0) k = std::stoi(argv[i+1]);
            else if (strcmp(argv[i], "-dk") == 0) deltaK = std::stoi(argv[i+1]);
            else if (strcmp(argv[i], "-ln") == 0) lNeg = std::stoi(argv[i+1]);
            else if (strcmp(argv[i], "-lp") == 0) lPoz = std::stoi(argv[i+1]);
            else if (strcmp(argv[i], "-o") == 0) outputFile = argv[i+1];
        }
        
        if (!generateInstance(n, k, deltaK, lNeg, lPoz, outputFile)) {
            std::cerr << "Failed to generate instance!\n";
            return 1;
        }
        std::cout << "Instance generated successfully: " << outputFile << std::endl;
        
    } else if (mode == "test_instance") {
        std::string inputFile;
        std::string outputFile = "results.txt";
        int processId = 0; // Domyślna wartość

        for (int i = 2; i < argc; i += 2) {
            if (i + 1 >= argc) break;
            
            if (strcmp(argv[i], "-i") == 0) inputFile = argv[i+1];
            else if (strcmp(argv[i], "-o") == 0) outputFile = argv[i+1];
            else if (strcmp(argv[i], "-pid") == 0) processId = std::stoi(argv[i+1]); // Nowy argument
        }
        
        if (inputFile.empty()) {
            std::cerr << "Input file must be specified for test_instance mode!\n";
            printUsage();
            return 1;
        }
        
        DNAInstance instance;
        if (!InstanceIO::loadInstance(inputFile, instance)) {
            std::cerr << "Failed to load instance from " << inputFile << std::endl;
            return 1;
        }
        
        runGeneticAlgorithm(instance, outputFile, processId); // Przekazanie processId
        
    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        printUsage();
        return 1;
    }

    return 0;
}
