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
DNAInstance generateInstance(int n, int k, int deltaK, int lNeg, int lPoz) {
    try {
        // Create both negative and positive error introducers
        auto negativeErrorIntroducer = ErrorIntroducerFactory::createNegativeErrorIntroducer(lNeg);
        auto positiveErrorIntroducer = ErrorIntroducerFactory::createPositiveErrorIntroducer(lPoz, k);

        DNAInstanceBuilder builder;
        builder.setN(n)
               .setK(k)
               .setDeltaK(deltaK)
               .setLNeg(lNeg)
               .setLPoz(lPoz)
               .setRepAllowed(true)
               .setProbablePositive(0)
               .buildDNA()
               .buildSpectrum()
               .applyError(negativeErrorIntroducer.get())  // Apply negative errors
               .applyError(positiveErrorIntroducer.get()); // Apply positive errors

        return builder.getInstance();
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to generate instance: " + std::string(e.what()));
        throw;
    }
}

// Declare the function (it's defined in genetic_algorithm_runner.cpp)
void runGeneticAlgorithm(const DNAInstance& instance,
                        const std::string& outputFile,
                        int processId,
                        const std::string& configFile,
                        bool debugMode);

double runGeneticAlgorithmWithFitness(const DNAInstance& instance,
                          [[maybe_unused]] const std::string& outputFile,
                          int maxIterations,
                          [[maybe_unused]] const std::string& logFile,
                          [[maybe_unused]] bool verbose) {
    GeneticConfig config;
    config.maxGenerations = maxIterations;
    config.populationSize = 100;
    config.mutationProbability = 0.1;
    config.crossoverProbability = 0.8;
    config.targetFitness = 0.95;
    config.tournamentSize = 5;

    auto representation = std::make_unique<DirectDNARepresentation>();
    GeneticAlgorithm ga(std::move(representation), config);

    std::string result = ga.run(instance);
    try {
        return std::stod(result);
    } catch (const std::exception& e) {
        std::cerr << "Error converting result to double: " << e.what() << std::endl;
        return 0.0;
    }
}

double evaluateParameterSet(const DNAInstance& testInstance, 
                        const ParameterSet& ps, 
                        [[maybe_unused]] double (*calculateFitness)(const DNAInstance&, const std::vector<int>&)) {
    GeneticConfig config;
    config.maxGenerations = 100;  // Default value
    config.populationSize = ps.contains("populationSize") ? ps.getInt("populationSize") : 100;
    config.mutationProbability = ps.contains("mutationRate") ? ps.getDouble("mutationRate") : 0.1;
    config.crossoverProbability = ps.contains("crossoverProb") ? ps.getDouble("crossoverProb") : 0.8;
    config.targetFitness = 0.95;  // Default value
    config.tournamentSize = 5;    // Default value

    auto representation = std::make_unique<DirectDNARepresentation>();
    GeneticAlgorithm ga(std::move(representation), config);

    std::string result = ga.run(testInstance);
    try {
        return std::stod(result);
    } catch (const std::exception& e) {
        std::cerr << "Error converting result to double: " << e.what() << std::endl;
        return 0.0;
    }
}

TuningResult runParameterTuning(const DNAInstance& testInstance, const std::vector<ParameterSet>& parameterSets) {
    ParameterSet bestSet;
    double bestFitness = std::numeric_limits<double>::lowest();

    for (const auto& ps : parameterSets) {
        double fitness = evaluateParameterSet(testInstance, ps, nullptr);  // Replace nullptr with actual fitness function
        if (fitness > bestFitness) {
            bestFitness = fitness;
            bestSet = ps;
        }
    }

    TuningResult tr;
    tr.parameterSet = bestSet;
    tr.fitness = bestFitness;
    tr.executionTime = 0.0;  // Add actual timing if needed
    return tr;
}

int main(int argc, char* argv[]) {
    // Initialize logger first
    Logger::initialize("log.txt");

    try {
        // Initialize variables with proper construction
        std::string mode;
        std::string configFile;
        bool debugMode = false;

        // Log program start with properly constructed strings
        LOG_INFO("Starting program with " + std::to_string(argc) + " arguments");

        // Log mode with properly constructed strings
        if (argc > 1) {
            mode = argv[1];  // Direct assignment is safe here
            LOG_INFO("Mode set to: " + mode);
        }

        // Log all arguments with properly constructed strings
        std::string argsMsg = "Arguments:";
        for (int i = 0; i < argc; ++i) {
            argsMsg += "\n  " + std::to_string(i) + ": " + argv[i];
        }
        LOG_INFO(argsMsg);

        // Parse arguments with proper string handling
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];  // Direct assignment
            if (arg == "-c" && i + 1 < argc) {
                configFile = argv[++i];  // Direct assignment
                LOG_INFO("Config file set to: " + configFile);
            } else if (arg == "-d") {
                debugMode = true;
                LOG_INFO("Debug mode enabled");
            }
        }

        // Validate config file path with proper error handling
        if (!configFile.empty()) {
            try {
                std::filesystem::path configPath(configFile);
                if (!std::filesystem::exists(configPath)) {
                    std::string errMsg = "Config file not found: " + configFile;
                    LOG_ERROR(errMsg);
                    return 1;
                }
                std::string existsMsg = "Config file exists: " + configFile;
                LOG_INFO(existsMsg);
            } catch (const std::exception& e) {
                std::string errMsg = "Invalid config file path: " + std::string(e.what());
                LOG_ERROR(errMsg);
                return 1;
            }
        }

        if (mode == "test_instance") {
            LOG_INFO("Entering test_instance mode");
            
            // Initialize all variables with proper construction
            std::string inputFile = "";
            std::string outputFile = "";
            int processId = -1;
            std::string difficulty = "";
            
            // Parse test_instance specific arguments with better error handling
            for (int i = 2; i < argc; ++i) {
                if (i + 1 >= argc) continue;  // Prevent buffer overrun
                
                const std::string arg = argv[i];
                const std::string value = argv[i + 1];
                
                if (arg == "-i") {
                    inputFile = value;
                    LOG_INFO("Input file set to: " + inputFile);
                    ++i;
                } else if (arg == "-o") {
                    outputFile = value;
                    LOG_INFO("Output file set to: " + outputFile);
                    ++i;
                } else if (arg == "-pid") {
                    try {
                        processId = std::stoi(value);
                        if (processId < 0) {
                            LOG_ERROR("Process ID must be non-negative");
                            return 1;
                        }
                        LOG_INFO("Process ID set to: " + std::to_string(processId));
                    } catch (const std::exception& e) {
                        LOG_ERROR("Invalid process ID: " + value + " - " + e.what());
                        return 1;
                    }
                    ++i;
                } else if (arg == "-diff") {
                    difficulty = value;
                    // Validate difficulty level
                    std::vector<std::string> validDifficulties = {"Easy", "Medium", "Hard"};
                    if (std::find(validDifficulties.begin(), validDifficulties.end(), difficulty) 
                        == validDifficulties.end()) {
                        LOG_ERROR("Invalid difficulty level: " + difficulty + ". Must be one of: Easy, Medium, Hard");
                        return 1;
                    }
                    LOG_INFO("Difficulty set to: " + difficulty);
                    ++i;
                }
            }
            
            // Validate required arguments
            std::vector<std::pair<std::string, std::string>> requiredArgs = {
                {"input file (-i)", inputFile},
                {"output file (-o)", outputFile},
                {"difficulty (-diff)", difficulty}
            };
            
            for (const auto& [name, value] : requiredArgs) {
                if (value.empty()) {
                    LOG_ERROR("Missing required argument: " + name);
                    return 1;
                }
            }
            
            if (processId == -1) {
                LOG_ERROR("Missing or invalid process ID (-pid)");
                return 1;
            }
            
            // Validate input file with proper path handling
            std::filesystem::path inputPath;
            try {
                inputPath = std::filesystem::path(inputFile);
                if (!std::filesystem::exists(inputPath)) {
                    LOG_ERROR("Input file not found: " + inputFile);
                    return 1;
                }
                LOG_INFO("Input file exists: " + inputFile);
            } catch (const std::exception& e) {
                LOG_ERROR("Invalid input file path: " + std::string(e.what()));
                return 1;
            }
            
            // Create output directory with proper error handling
            std::filesystem::path outputPath;
            try {
                outputPath = std::filesystem::path(outputFile);
                auto outputDir = outputPath.parent_path();
                if (!outputDir.empty()) {
                    std::error_code ec;
                    if (!std::filesystem::exists(outputDir)) {
                        if (!std::filesystem::create_directories(outputDir, ec)) {
                            LOG_ERROR("Failed to create output directory: " + 
                                     outputDir.string() + " - " + ec.message());
                            return 1;
                        }
                        LOG_INFO("Created output directory: " + outputDir.string());
                    }
                }
            } catch (const std::exception& e) {
                LOG_ERROR("Invalid output file path: " + std::string(e.what()));
                return 1;
            }
            
            // Load and validate instance with proper error handling
            DNAInstance instance;
            try {
                LOG_INFO("Loading instance from: " + inputFile);
                // Load instance from file
                try {
                    instance = InstanceIO::loadInstance(inputFile);
                } catch (const std::exception& e) {
                    std::cerr << "Error loading instance: " << e.what() << std::endl;
                    return 1;
                }
                
                // Validate instance
                if (instance.getSpectrum().empty() || instance.getDNA().empty() || instance.getK() <= 0) {
                    LOG_ERROR("Invalid instance: missing required data");
                    return 1;
                }
                
                LOG_INFO("Instance loaded and validated successfully");
                LOG_INFO("Instance details: N=" + std::to_string(instance.getN()) + 
                        ", K=" + std::to_string(instance.getK()));
            } catch (const std::exception& e) {
                LOG_ERROR("Exception while loading instance: " + std::string(e.what()));
                return 1;
            }
            
            // Run genetic algorithm with proper error handling
            try {
                LOG_INFO("Starting genetic algorithm with process ID: " + std::to_string(processId));
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

            if (!InstanceIO::saveInstance(generateInstance(n, k, deltaK, lNeg, lPoz), outputFile)) {
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

            // Add helper function for fitness calculation
            auto calculateFitness = [](const std::string& solution, const DNAInstance& instance) -> double {
                // Calculate Hamming distance between solution and target DNA
                const std::string& targetDNA = instance.getDNA();
                if (solution.length() != targetDNA.length()) {
                    return 0.0;  // Invalid solution
                }
                
                int matches = 0;
                for (size_t i = 0; i < solution.length(); ++i) {
                    if (solution[i] == targetDNA[i]) {
                        matches++;
                    }
                }
                
                return static_cast<double>(matches) / solution.length();
            };
            
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
            auto evaluateFunc = [&testInstance, &calculateFitness](const ParameterSet& ps) {
                // Create local config with default values
                GeneticConfig config;
                
                // Update config with parameters from parameter set
                if (ps.contains("populationSize")) {
                    config.populationSize = ps.getInt("populationSize");
                }
                if (ps.contains("mutationRate")) {
                    config.mutationProbability = ps.getDouble("mutationRate");
                }
                if (ps.contains("crossoverProb")) {
                    config.crossoverProbability = ps.getDouble("crossoverProb");
                }
                
                // Create genetic algorithm
                auto representation = std::make_unique<DirectDNARepresentation>();
                GeneticAlgorithm ga(std::move(representation), config);
                
                // Run the algorithm
                auto startTime = std::chrono::high_resolution_clock::now();
                auto solution = ga.run(testInstance);
                auto endTime = std::chrono::high_resolution_clock::now();
                
                // Calculate fitness and create result
                TuningResult tr;
                tr.parameterSet = ps;
                tr.fitness = calculateFitness(solution, testInstance);
                tr.executionTime = std::chrono::duration<double>(endTime - startTime).count();
                tr.extraMetrics["converged"] = 1.0;
                return tr;
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

                // Create a new config instance with default values
                GeneticConfig config;

                // Update parameters
                for (const auto &[key, value] : ps.params) {
                    if (key == "populationSize") {
                        config.populationSize = std::stoi(value);
                    } else if (key == "mutationRate") {
                        config.mutationProbability = std::stod(value);
                    }
                }

                // Run hybrid tuning
                HybridOnePlusLambdaEA hybrid;
                ParameterSet bestSet = hybrid.runHybridOnePlusLambdaEA(ps, instance);

                TuningResult tr;
                tr.parameterSet = bestSet;
                tr.fitness = config.targetFitness;
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

        // Cleanup logger before exit
        Logger::cleanup();
        return 0;
    } catch (const std::exception& e) {
        LOG_ERROR("Unhandled exception: " + std::string(e.what()));
        Logger::cleanup();
        return 1;
    } catch (...) {
        LOG_ERROR("Unknown error occurred");
        Logger::cleanup();
        return 1;
    }
}
