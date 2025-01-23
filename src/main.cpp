#include <iostream>
#include <random>
#include <chrono>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>

// Pliki projektu
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

// Nowy include do konfiguracji AG:
#include "configuration/genetic_algorithm_configuration.h"

#include "metaheuristics/adaptive_crossover.h"
#include "benchmark/adaptive_crossover_benchmark.h"

// Funkcja do obliczania odległości Levenshteina (może być również w utils)
int levenshteinDistance(const std::string &s1, const std::string &s2) {
    int len1 = static_cast<int>(s1.size());
    int len2 = static_cast<int>(s2.size());
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));

    for(int i = 0; i <= len1; ++i) dp[i][0] = i;
    for(int j = 0; j <= len2; ++j) dp[0][j] = j;

    for(int i = 1; i <= len1; ++i) {
        for(int j = 1; j <= len2; ++j) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            dp[i][j] = std::min({
                dp[i-1][j] + 1,
                dp[i][j-1] + 1,
                dp[i-1][j-1] + cost
            });
        }
    }
    return dp[len1][len2];
}

// Funkcja do wypisywania użycia programu:
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
              << "  -pid <value>      - Process ID (unique identifier)\n"
              << "  -diff <string>    - Difficulty/test name (e.g. Easy, Medium, Hard)\n\n"
              << "Additional option for any mode:\n"
              << "  -cfg <file>       - Path to GA config file (key=value style)\n";
}

// Funkcja uruchamiająca Algorytm Genetyczny
void runGeneticAlgorithm(const DNAInstance& instance,
                         const std::string& outputFile = "",
                         int processId = 0,
                         const std::string& difficulty = "Unknown")
{
    // Wyciągamy singleton GAConfig (już wczytany z pliku przez main)
    auto& config = GAConfig::getInstance();

    // Ustawiamy cache, jeśli chcemy
    auto cache = std::make_shared<CachedPopulation>();
    config.setCache(cache);

    // Tworzymy obiekt GA
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        config.getCrossover(config.crossoverType),
        config.getMutation(),
        config.getReplacement(),
        config.getFitness(),
        config.getStopping(),
        cache
    );

    // Ustawienie callbacku do aktualizacji postępu (opcjonalne)
    ga.setProgressCallback([processId](int generation,
                                       int maxGenerations,
                                       double bestFitness,
                                       double coverage,
                                       double edgeScore,
                                       double theoreticalMax)
    {
        double progress = (static_cast<double>(generation) / maxGenerations) * 100.0;
        std::ostringstream status;
        status << "Gen " << generation << "/" << maxGenerations
               << " Fit: " << bestFitness
               << " Cov: " << coverage
               << " Edge: " << edgeScore;

        // Format: PROGRESS_UPDATE:PID:progress:status:bestFitness:coverage:edgeScore:theoreticalMax
        std::cout << "PROGRESS_UPDATE:" << processId << ":"
                  << std::fixed << std::setprecision(2) << progress << ":"
                  << status.str() << ":"
                  << bestFitness << ":"
                  << coverage << ":"
                  << edgeScore << ":"
                  << theoreticalMax
                  << std::endl;
        std::cout.flush();
    });

    // Identyfikator procesu
    ga.setProcessId(processId);

    // Uruchamiamy GA
    ga.run(instance);

    // Odczytujemy najlepsze znalezione rozwiązanie
    std::string reconstructedDNA = ga.getBestDNA();
    std::string originalDNA = instance.getDNA();

    int distance = levenshteinDistance(originalDNA, reconstructedDNA);

    // Zapis do pliku lub wypis w stdout
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
            std::cout << "PROGRESS_UPDATE:" << processId << ":100:Completed:0:0:0:0\n";
            std::cout.flush();
        } else {
            std::cerr << "Failed to open output file: " << outputFile << std::endl;
        }
    }
}

// Funkcja pomocnicza do generowania instancji
bool generateInstance(int n,
                      int k,
                      int deltaK,
                      int lNeg,
                      int lPoz,
                      const std::string& outputFile)
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
    // Dodajemy prostą obsługę parametru -cfg:
    std::string configFile = "config.cfg"; // domyślna ścieżka do pliku z parametrami GA

    // Parsowanie argumentów w trybie -cfg
    for (int i = 2; i < argc; ++i) {
        if (std::strcmp(argv[i], "-cfg") == 0) {
            if (i + 1 < argc) {
                configFile = argv[i + 1];
                ++i; // pomijamy wartość
            } else {
                std::cerr << "Error: -cfg requires a file path.\n";
                return 1;
            }
        }
    }

    // Wczytujemy parametry GA z pliku, np. "config.cfg"
    auto& config = GAConfig::getInstance();
    if (!config.loadFromFile(configFile)) {
        std::cerr << "Failed to load GA configuration from " << configFile << "\n";
        return 1;
    }

    // Teraz przechodzimy do logiki poszczególnych trybów:
    if (mode == "debug") {
        // Domyślne parametry do generowania instancji
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

        // Uruchamiamy GA na tej instancji
        runGeneticAlgorithm(instance, "");
        // Jeśli używasz profilera, upewnij się, że jest on poprawnie zaimplementowany
        // Profiler::getInstance().saveReport("debug_profiling_report.csv");

    } else if (mode == "generate_instance") {
        // Parametry domyślne
        int n = 400, k = 8, deltaK = 1, lNeg = 10, lPoz = 10;
        std::string outputFile = "generated_instance.txt";

        // Parsowanie argumentów w parze: -n, -k, -dk, -ln, -lp, -o
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
                // Już obsłużone wcześniej
                ++i;
            }
        }

        if (!generateInstance(n, k, deltaK, lNeg, lPoz, outputFile)) {
            std::cerr << "Failed to generate instance!\n";
            return 1;
        }
        std::cout << "Instance generated successfully: " << outputFile << std::endl;

    } else if (mode == "test_instance") {
        std::string inputFile;
        std::string outputFile = "results.txt";
        int processId = 0;
        std::string difficulty = "Unknown";

        // Parsowanie argumentów w parze: -i, -o, -pid, -diff
        for (int i = 2; i < argc; ++i) {
            if (std::strcmp(argv[i], "-i") == 0) {
                if (i + 1 < argc) {
                    inputFile = argv[i + 1];
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-o") == 0) {
                if (i + 1 < argc) {
                    outputFile = argv[i + 1];
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-pid") == 0) {
                if (i + 1 < argc) {
                    processId = std::stoi(argv[i + 1]);
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-diff") == 0) {
                if (i + 1 < argc) {
                    difficulty = argv[i + 1];
                    ++i;
                }
            } else if (std::strcmp(argv[i], "-cfg") == 0) {
                // Już obsłużone wcześniej
                ++i;
            }
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

        std::cout << "CURRENT_TEST:"
                  << processId << ":"
                  << difficulty << ":"
                  << inputFile
                  << std::endl;

        // Uruchomienie GA na wczytanej instancji
        runGeneticAlgorithm(instance, outputFile, processId, difficulty);

    } else {
        std::cerr << "Unknown mode: " << mode << std::endl;
        printUsage();
        return 1;
    }

    return 0;
}
