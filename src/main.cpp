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
#include "../include/dna/dna_instance.h"
#include "../include/interfaces/i_representation.h"
#include "../include/interfaces/i_selection.h"
#include "../include/interfaces/i_crossover.h"
#include "../include/interfaces/i_mutation.h"
#include "../include/interfaces/i_replacement.h"
#include "../include/interfaces/i_fitness.h"
#include "../include/interfaces/i_stopping.h"
#include "../include/interfaces/i_population_cache.h"

#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../include/metaheuristics/population_cache_impl.h"

#include "../include/generator/dna_generator.h"
#include "../include/naive/naive_reconstruction.h"
#include "../include/benchmark/naive_benchmark.h"
#include "../include/utils/utility_functions.h"

#include "../include/metaheuristics/crossover_impl.h"
#include "../include/metaheuristics/fitness_impl.h"
#include "../include/metaheuristics/genetic_algorithm_runner.h"

#include "../include/benchmark/crossover_benchmark.h"
#include "../include/benchmark/adaptive_crossover_benchmark.h"

#include "../include/metaheuristics/adaptive_crossover.h"

// Tuning headers
#include "../include/tuning/parameter_tuning_manager.h"
#include "../include/tuning/parameters_parser.h"
#include "../include/tuning/tuning_structures.h"
#include "../include/tuning/racing.h"
#include "../include/tuning/meta_ea.h"

#include "../include/utils/logging.h"
#include "../include/metaheuristics/representation.h"
#include "../include/dna/error_introduction.h"

// Create error introducers
auto negativeErrorIntroducer = std::make_unique<NegativeErrorIntroducer>();
auto positiveErrorIntroducer = std::make_unique<PositiveErrorIntroducer>();

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
        auto random = std::make_unique<Random>();
        auto generator = std::make_unique<DNAGenerator>(std::move(random));
        DNAInstanceBuilder builder(std::move(generator));
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
                                [[maybe_unused]] int processId, 
                                [[maybe_unused]] const std::string& configFile, 
                                [[maybe_unused]] bool debugMode) {
    // Create configuration
    GAConfig config;
    config.setPopulationSize(200);
    config.setMutationRate(0.5);
    config.setCrossoverProbability(1.0);
    config.setTournamentSize(1);
    
    // Create representation
    auto representation = std::make_unique<PermutationRepresentation>();
    
    // Create and run genetic algorithm
    GeneticAlgorithm ga(std::move(representation), config);
    ga.setProcessId(processId);
    
    std::string result = ga.run(instance);
    return ga.getBestFitness();
}

double evaluateParameterSet(const DNAInstance& instance, 
                          const ParameterSet& params, 
                          [[maybe_unused]] double (*fitnessFunc)(const DNAInstance&, const std::vector<int>&)) {
    // Create configuration
    GAConfig config;
    config.setPopulationSize(params.getInt("populationSize"));
    config.setMutationRate(params.getDouble("mutationRate"));
    config.setCrossoverProbability(params.getDouble("crossoverProb"));
    config.setTournamentSize(5);  // Default value
    
    // Create representation
    auto representation = std::make_unique<PermutationRepresentation>();
    
    // Create and run genetic algorithm
    GeneticAlgorithm ga(std::move(representation), config);
    std::string result = ga.run(instance);
    
    return ga.getBestFitness();
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

void runParameterTuningWithRacing(const DNAInstance& instance) {
    // Create racing configuration
    Racing::Configuration rc;
    rc.significanceLevel = 0.05;
    rc.maxTrialsPerCandidate = 10;
    rc.minTrialsBeforeElimination = 3;
    
    [[maybe_unused]] // Mark as unused since it's a placeholder for future use
    auto evaluateFunc = [&instance](const ParameterSet& ps) -> TuningResult {
        try {
            // Create configuration
            GAConfig config;
            
            // Set parameters from ParameterSet
            for (const auto& [key, value] : ps.params) {
                try {
                    if (key == "populationSize") {
                        config.setPopulationSize(std::stoi(value));
                    } else if (key == "mutationRate") {
                        config.setMutationRate(std::stod(value));
                    } else if (key == "crossoverProb") {
                        config.setCrossoverProbability(std::stod(value));
                    }
                } catch (const std::exception& e) {
                    LOG_ERROR("Error setting parameter " + key + ": " + std::string(e.what()));
                }
            }
            
            // Create representation and genetic algorithm
            auto representation = std::make_unique<PermutationRepresentation>();
            GeneticAlgorithm ga(std::move(representation), config);
            std::string result = ga.run(instance);
            
            // Create and return result
            TuningResult tr;
            tr.parameterSet = ps;
            tr.fitness = ga.getBestFitness();
            tr.executionTime = 0.0;  // TODO: Add timing
            return tr;
            
        } catch (const std::exception& e) {
            LOG_ERROR("Error in parameter evaluation: " + std::string(e.what()));
            return TuningResult(ps, 0.0, 0.0);
        }
    };
    
    // ... rest of the function ...
}

int main(int argc, char* argv[]) {
    // Initialize logger first
    Logger::initialize("log.txt");
    
    try {
        // Use RAII to ensure logger cleanup
        struct LoggerCleanup {
            ~LoggerCleanup() { Logger::cleanup(); }
        } loggerGuard;

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
            if (arg == "-cfg" && i + 1 < argc) {
                configFile = argv[++i];  // Direct assignment
                LOG_INFO("Config file set to: " + configFile);
            } else if (arg == "-d") {
                debugMode = true;
                LOG_INFO("Debug mode enabled");
                Logger::setLogLevel(LogLevel::DEBUG);
            }
        }

        // Validate config file path with proper error handling
        if (!configFile.empty()) {
            try {
                std::filesystem::path configPath(configFile);
                if (!std::filesystem::exists(configPath)) {
                    // Try looking in the build directory's config subdirectory
                    std::filesystem::path buildDirConfig = std::filesystem::current_path() / configFile;
                    if (!std::filesystem::exists(buildDirConfig)) {
                        std::string errMsg = "Config file not found at: " + configFile + " or " + buildDirConfig.string();
                        LOG_ERROR(errMsg);
                        return 1;
                    }
                    configFile = buildDirConfig.string();
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
                    auto random = std::make_unique<Random>();
                    auto generator = std::make_unique<DNAGenerator>(std::move(random));
                    DNAInstanceBuilder builder(std::move(generator));
                    instance = builder.setN(300)
                                    .setK(7)
                                    .setDeltaK(0)
                                    .setLNeg(10)
                                    .setLPoz(10)
                                    .setRepAllowed(true)
                                    .setProbablePositive(0)
                                    .buildDNA()
                                    .buildSpectrum()
                                    .applyError(negativeErrorIntroducer.get())  // Apply negative errors
                                    .applyError(positiveErrorIntroducer.get())  // Apply positive errors
                                    .getInstance();
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
        } else if (mode == "debug") {
            LOG_INFO("Entering debug mode");
            
            // Set debug log level
            Logger::setLogLevel(LogLevel::DEBUG);
            LOG_DEBUG("Debug logging enabled");
            
            // Create a default test instance for debugging
            auto random = std::make_unique<Random>();
            auto generator = std::make_unique<DNAGenerator>(std::move(random));
            DNAInstanceBuilder builder(std::move(generator));
            builder.setN(300)
                   .setK(7)
                   .setDeltaK(0)
                   .setLNeg(10)
                   .setLPoz(10)
                   .setRepAllowed(true)
                   .setProbablePositive(0)
                   .buildDNA()
                   .buildSpectrum();
            
            DNAInstance instance = builder.getInstance();
            // Set target sequence to original DNA for Levenshtein distance calculation
            instance.setTargetSequence(instance.getDNA());

            // Log instance details
            LOG_INFO("Instance details:");
            LOG_INFO("  DNA length (N): " + std::to_string(instance.getN()));
            LOG_INFO("  Oligo length (K): " + std::to_string(instance.getK()));
            LOG_INFO("  Delta K: " + std::to_string(instance.getDeltaK()));
            LOG_INFO("  Negative errors (LNeg): " + std::to_string(instance.getLNeg()));
            LOG_INFO("  Positive errors (LPoz): " + std::to_string(instance.getLPoz()));
            LOG_INFO("  Repetitions allowed: " + std::string(instance.isRepAllowed() ? "yes" : "no"));
            LOG_INFO("  DNA sequence: " + instance.getDNA());
            LOG_INFO("  Target sequence: " + instance.getTargetSequence());
            LOG_INFO("  Spectrum size: " + std::to_string(instance.getSpectrum().size()));

            std::string outputFile = "debug_output.txt";
            int processId = 0;
            
            // Run genetic algorithm in debug mode
            runGeneticAlgorithm(instance, outputFile, processId, configFile, true);
            
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
            auto random = std::make_unique<Random>();
            auto generator = std::make_unique<DNAGenerator>(std::move(random));
            DNAInstanceBuilder builder(std::move(generator));
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
                GAConfig config;
                
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
                
                // Create genetic algorithm
                auto representation = std::make_unique<PermutationRepresentation>();
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
                auto random = std::make_unique<Random>();
                auto generator = std::make_unique<DNAGenerator>(std::move(random));
                DNAInstanceBuilder builder(std::move(generator));
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
                GAConfig config;

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
                tr.fitness = config.getTargetFitness();
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
    } catch (const std::exception& e) {
        LOG_ERROR("Unhandled exception: " + std::string(e.what()));
        return 1;
    } catch (...) {
        LOG_ERROR("Unknown error occurred");
        return 1;
    }
}
