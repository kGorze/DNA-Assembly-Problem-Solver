#include <iostream>
#include <random>
#include <chrono>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>
#include <ctime>
#include <vector>
#include <thread>
#include <mutex>
#include <filesystem>
#include <algorithm>  // for std::transform

// Project headers
#include "dna/dna_instance.h"
#include "interfaces/i_representation.h"
#include "interfaces/i_selection.h"
#include "interfaces/i_crossover.h"
#include "interfaces/i_mutation.h"
#include "interfaces/i_replacement.h"
#include "interfaces/i_fitness.h"
#include "interfaces/i_stopping.h"
#include "interfaces/i_population_cache.h"

#include "configuration/genetic_algorithm_configuration.h"
#include "metaheuristics/population_cache_impl.h"

#include "generator/dna_generator.h"
#include "naive/naive_reconstruction.h"
#include "benchmark/naive_benchmark.h"
#include "utils/utility_functions.h"

#include "metaheuristics/crossover_impl.h"
#include "metaheuristics/fitness_impl.h"
#include "metaheuristics/genetic_algorithm_runner.h"

#include "benchmark/crossover_benchmark.h"
#include "benchmark/adaptive_crossover_benchmark.h"

#include "metaheuristics/adaptive_crossover.h"

// Tuning headers
#include "tuning/parameter_tuning_manager.h"
#include "tuning/parameters_parser.h"
#include "tuning/tuning_structures.h"
#include "tuning/racing.h"
#include "tuning/meta_ea.h"

#include "utils/logging.h"


// Funkcja do wypisywania użycia programu:
void printUsage() {
    std::cout << "Usage: dna_reconstruction <mode> [options]\n\n"
              << "Modes:\n"
              << "  debug             - Run in debug mode with default settings\n"
              << "  generate_instance - Generate test instances\n"
              << "  test_instance     - Solve DNA reconstruction from input file\n\n"
              << "  tuning            - Run parameter tuning\n"
              << "  tuning_hybrid     - Run advanced hybrid (1+lambda) ES + Racing\n\n"
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
              << "  -pid <value>      - Process ID (unique identifier)\n"
              << "  -diff <string>    - Difficulty/test name (e.g. Easy, Medium, Hard)\n\n"
              << "Options for tuning:\n"
              << "  -cfg <file>       - Path to GA config file (key=value style)\n"
              << "  -out <file>       - Output CSV file for tuning results\n\n"
              << "Options for tuning_hybrid:\n"
              << "  (same as tuning, but uses the new HybridOnePlusLambdaEA approach)\n\n"
              << "Additional option for any mode:\n"
              << "  -cfg <file>       - Path to GA config file (key=value style)\n";
}

// Funkcja pomocnicza do generowania instancji
bool generateInstance(int n, int k, int deltaK, int lNeg, int lPoz, const std::string &outputFile)
{
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

    DNAInstance instance = builder.getInstance();
    
    // Introduce errors AFTER setting start index
    auto startFrag = instance.getDNA().substr(0, k);
    const auto& spectrum = instance.getSpectrum();
    
    int startIdx = -1;
    for (int i = 0; i < (int)spectrum.size(); i++) {
        if (spectrum[i] == startFrag) {
            startIdx = i;
            break;
        }
    }
    
    instance.setStartIndex(startIdx);

    // Now introduce errors
    if (lNeg > 0) {
        NegativeErrorIntroducer negErr(lNeg);
        negErr.introduceErrors(instance);
    }

    if (lPoz > 0) {
        PositiveErrorIntroducer posErr(lPoz);
        posErr.introduceErrors(instance);
    }

    return InstanceIO::saveInstance(instance, outputFile);
}

// Declare the function (it's defined in genetic_algorithm_runner.cpp)
void runGeneticAlgorithm(const DNAInstance& instance,
                        const std::string& outputFile,
                        int processId,
                        const std::string& configFile,
                        bool debugMode);

int main(int argc, char* argv[]) {
    LOG_INFO("Starting program with " + std::to_string(argc) + " arguments");
    
    if (argc < 2) {
        LOG_ERROR("No mode specified (argc < 2)");
        printUsage();
        return 1;
    }

    std::string mode = argv[1];
    LOG_INFO("Mode: " + mode);
    
    // Domyślna ścieżka do pliku config
    std::string configFile = "config.cfg";
    bool debugMode = false;

    // Log all arguments
    for (int i = 0; i < argc; ++i) {
        LOG_INFO("argv[" + std::to_string(i) + "]: " + std::string(argv[i]));
    }

    // Parsujemy argumenty w poszukiwaniu -cfg i --debug
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "-cfg") == 0) {
            if (i + 1 < argc) {
                configFile = argv[i + 1];
                LOG_INFO("Config file path set to: " + configFile);
                ++i; // pomijamy wartość
            } else {
                LOG_ERROR("Error: -cfg requires a file path");
                return 1;
            }
        } else if (std::strcmp(argv[i], "--debug") == 0) {
            debugMode = true;
            LOG_INFO("Debug output enabled");
        }
    }

    // Validate config file
    if (!std::filesystem::exists(configFile)) {
        LOG_ERROR("Config file not found: " + configFile);
        return 1;
    }
    LOG_INFO("Config file exists: " + configFile);

    if (mode == "test_instance") {
        LOG_INFO("Entering test_instance mode");
        
        std::string inputFile;
        std::string outputFile;
        int processId = -1;
        std::string difficulty;
        
        // Parse test_instance specific arguments
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
                inputFile = argv[i + 1];
                LOG_INFO("Input file set to: " + inputFile);
                ++i;
            } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
                outputFile = argv[i + 1];
                LOG_INFO("Output file set to: " + outputFile);
                ++i;
            } else if (std::strcmp(argv[i], "-pid") == 0 && i + 1 < argc) {
                try {
                    processId = std::stoi(argv[i + 1]);
                    LOG_INFO("Process ID set to: " + std::to_string(processId));
                } catch (const std::exception& e) {
                    LOG_ERROR("Invalid process ID: " + std::string(argv[i + 1]));
                    return 1;
                }
                ++i;
            } else if (std::strcmp(argv[i], "-diff") == 0 && i + 1 < argc) {
                difficulty = argv[i + 1];
                LOG_INFO("Difficulty set to: " + difficulty);
                ++i;
            }
        }
        
        // Validate required arguments
        if (inputFile.empty()) {
            LOG_ERROR("No input file specified (-i)");
            return 1;
        }
        if (outputFile.empty()) {
            LOG_ERROR("No output file specified (-o)");
            return 1;
        }
        if (processId == -1) {
            LOG_ERROR("No process ID specified (-pid)");
            return 1;
        }
        if (difficulty.empty()) {
            LOG_ERROR("No difficulty specified (-diff)");
            return 1;
        }
        
        // Validate input file exists
        if (!std::filesystem::exists(inputFile)) {
            LOG_ERROR("Input file not found: " + inputFile);
            return 1;
        }
        LOG_INFO("Input file exists: " + inputFile);
        
        // Create output directory if needed
        auto outputDir = std::filesystem::path(outputFile).parent_path();
        if (!outputDir.empty()) {
            std::error_code ec;
            std::filesystem::create_directories(outputDir, ec);
            if (ec) {
                LOG_ERROR("Failed to create output directory: " + outputDir.string() + " - " + ec.message());
                return 1;
            }
            LOG_INFO("Created output directory: " + outputDir.string());
        }
        
        // Load and validate instance
        DNAInstance instance;
        try {
            LOG_INFO("Loading instance from: " + inputFile);
            if (!InstanceIO::loadInstance(inputFile, instance)) {
                LOG_ERROR("Failed to load instance from file");
                return 1;
            }
            LOG_INFO("Instance loaded successfully");
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to load instance: " + std::string(e.what()));
            return 1;
        }
        
        // Run genetic algorithm
        try {
            LOG_INFO("Starting genetic algorithm");
            runGeneticAlgorithm(instance, outputFile, processId, configFile, debugMode);
            LOG_INFO("Genetic algorithm completed successfully");
            return 0;
        } catch (const std::exception& e) {
            LOG_ERROR("Genetic algorithm failed: " + std::string(e.what()));
            return 1;
        }
    } else if (mode == "generate_instance") {
        int n = 400, k = 8, deltaK = 1, lNeg = 10, lPoz = 10;
        std::string outputFile = "generated_instance.txt";

        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "-n") == 0) {
                if (i + 1 < argc) {
                    n = std::stoi(argv[i + 1]);
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-k") == 0) {
                if (i + 1 < argc) {
                    k = std::stoi(argv[i + 1]);
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-dk") == 0) {
                if (i + 1 < argc) {
                    deltaK = std::stoi(argv[i + 1]);
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-ln") == 0) {
                if (i + 1 < argc) {
                    lNeg = std::stoi(argv[i + 1]);
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-lp") == 0) {
                if (i + 1 < argc) {
                    lPoz = std::stoi(argv[i + 1]);
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-o") == 0) {
                if (i + 1 < argc) {
                    outputFile = argv[i + 1];
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-cfg") == 0) {
                ++i; // skip
            }
        }

        if (!generateInstance(n, k, deltaK, lNeg, lPoz, outputFile)) {
            std::cerr << "Failed to generate instance!\n";
            return 1;
        }
        std::cout << "Instance generated successfully: " << outputFile << std::endl;

    } else if (mode == "tuning") {
        std::string tuningOutputFile = "tuning_results.csv";
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "-out") == 0 && i + 1 < argc) {
                tuningOutputFile = argv[i + 1];
                ++i;
            }
        }

        // Create parameter tuning manager with output file
        ParameterTuningManager tuner(tuningOutputFile);
        
        // Generate candidate parameter sets
        std::vector<ParameterSet> candidateParams;
        
        // Add some example parameter sets
        ParameterSet ps1;
        ps1.params["populationSize"] = "100";
        ps1.params["mutationRate"] = "0.1";
        ps1.params["crossoverProb"] = "0.8";
        candidateParams.push_back(ps1);
        
        ParameterSet ps2;
        ps2.params["populationSize"] = "200";
        ps2.params["mutationRate"] = "0.2";
        ps2.params["crossoverProb"] = "0.7";
        candidateParams.push_back(ps2);
        
        // Set up racing configuration
        Racing::Configuration rc;
        rc.significanceLevel = 0.05;
        rc.maxTrialsPerCandidate = 10;
        rc.minTrialsBeforeElimination = 3;
        
        // Create test instance for evaluation
        DNAInstanceBuilder builder;
        builder.setN(300)
               .setK(7)
               .setDeltaK(2)
               .setLNeg(0)
               .setLPoz(0)
               .setRepAllowed(true)
               .buildDNA()
               .buildSpectrum();
        DNAInstance testInstance = builder.getInstance();
        
        // Define evaluation function
        auto evaluateFunc = [&testInstance](const ParameterSet& ps) {
            // Create local config for this evaluation
            GAConfig config;
            if (!config.loadFromFile("config.cfg")) {
                throw std::runtime_error("Failed to load configuration");
            }
            
            // Update config with parameters from parameter set
            if (ps.contains("populationSize")) {
                config.setPopulationSize(ps.getInt("populationSize"));
            }
            if (ps.contains("mutationRate")) {
                config.setMutationRate(ps.getDouble("mutationRate"));
            }
            if (ps.contains("crossoverProb")) {
                config.setCrossoverProbability(ps.getDouble("crossoverProb"));
            }
            
            // Create and run GA
            GeneticAlgorithm ga(
                config.getRepresentation(),
                config.getSelection(),
                config.getCrossover("order"), // Specify crossover type
                config.getMutation(),
                config.getReplacement(),
                config.getFitness(),
                config.getStopping(),
                config.getCache(),
                config
            );
            
            ga.run(testInstance);
            
            TuningResult result;
            result.parameterSet = ps;
            result.fitness = ga.getBestFitness();
            return result;
        };
        
        // Run racing
        tuner.runRacingOnly(candidateParams, rc, evaluateFunc);
        
        std::cout << "Parameter tuning completed. Results saved to " << tuningOutputFile << std::endl;
        return 0;

    } else if (mode == "tuning_hybrid") {
        std::string tuningOutputFile = "tuning_hybrid_results.csv";
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "-out") == 0) {
                if (i + 1 < argc) {
                    tuningOutputFile = argv[i + 1];
                    ++i;
                } else {
                    std::cerr << "Error: -out requires a filename.\n";
                    return 1;
                }
            }
        }

        ParameterTuningManager tuner(tuningOutputFile);

        Racing::Configuration rc;
        rc.significanceLevel = 0.05;
        rc.maxTrialsPerCandidate = 10;
        rc.minTrialsBeforeElimination = 3;
        rc.useBootstrap = true;
        rc.bootstrapSamples = 1000;

        auto evaluateFunc = [&](const ParameterSet &ps) -> TuningResult {
            // Create test instance
            DNAInstanceBuilder builder;
            builder.setN(300)
                   .setK(7)
                   .setDeltaK(2)
                   .setLNeg(0)
                   .setLPoz(0)
                   .setRepAllowed(true)
                   .buildDNA()
                   .buildSpectrum();
            DNAInstance instance = builder.getInstance();

            // Create a new config instance
            GAConfig config;
            if (!config.loadFromFile("config.cfg")) {
                throw std::runtime_error("Failed to load configuration");
            }

            // Update parameters
            for (const auto &[key, value] : ps.params) {
                if (key == "populationSize") {
                    config.setPopulationSize(std::stoi(value));
                } else if (key == "mutationRate") {
                    config.setMutationRate(std::stod(value));
                }
            }

            // Run hybrid tuning
            HybridOnePlusLambdaEA hybrid;
            ParameterSet bestSet = hybrid.runHybridOnePlusLambdaEA(ps, instance);

            TuningResult tr;
            tr.parameterSet = bestSet;
            tr.fitness = config.getGlobalBestFitness();
            tr.executionTime = 0.0;
            return tr;
        };

        // Generate candidates
        std::vector<ParameterSet> candidateParams = ParameterParser::generateGridOfCandidatesWithout();

        // Run racing
        tuner.runRacingOnly(candidateParams, rc, evaluateFunc);

        std::cout << "Hybrid (1+lambda) tuning completed. Results: " << tuningOutputFile << "\n";
    }
    else {
        LOG_ERROR("Unknown mode: " + mode);
        printUsage();
        return 1;
    }

    return 0;
}
