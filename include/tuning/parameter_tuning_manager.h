//
// Created by konrad_guest on 23/01/2025.
//

#ifndef PARAMETER_TUNING_MANAGER_H
#define PARAMETER_TUNING_MANAGER_H

// System includes
#include <fstream>
#include <functional>
#include <vector>
#include <string>
#include <chrono>

// Project includes
#include "racing.h"
#include "parameters_parser.h"
#include "tuning_structures.h"
#include "meta_ea.h"
#include "metaheuristics/genetic_algorithm.h"
#include "../generator/dna_generator.h"
#include "dna/dna_instance_builder.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "utils/logging.h"
#include <random>
#include <algorithm>
#include "../dna/dna_instance.h"

// Forward declaration of runGeneticAlgorithm
void runGeneticAlgorithm(const DNAInstance& instance,
                        const std::string& outputFile,
                        int processId,
                        const std::string& configFile,
                        bool debugMode);

/**
 * Główna klasa spinająca cały proces tuningu.
 * Można jej użyć do:
 *   - Samodzielnego trybu "tylko Racing" (podajemy listę kandydatów).
 *   - Trybu "Meta-EA + Racing" (szukanie optymalnych parametrów w dużej przestrzeni).
 */
class ParameterTuningManager {
public:
    explicit ParameterTuningManager(const std::string& outputFile) : m_outputFile(outputFile) {}

    /**
     * Przykładowy interfejs do "tylko Racing".
     * @param candidateParameters - lista kandydatów
     * @param racingCfg - parametry Racing
     */
    void runRacingOnly(const std::vector<ParameterSet> &candidateParameters,
                       const Racing::Configuration &racingCfg)
    {
        Racing::Manager rm(racingCfg);

        // Funkcja-lambda wywołująca Twój AE z danym zestawem parametrów.
        auto evaluateFunc = [this](const ParameterSet &ps) {
            return this->runOneEvaluation(ps);
        };

        auto results = rm.runRacing(candidateParameters, evaluateFunc);

        // Zapisz do CSV
        saveResults(results);
    }

    std::vector<TuningResult> runRacingOnly(
        const std::vector<ParameterSet>& candidateParameters,
        const Racing::Configuration& config,
        std::function<TuningResult(const ParameterSet&)> evaluateFunc);

    /**
     * Tryb "Meta-EA + Racing".
     * @param metaCfg - parametry Meta-EA (np. populacja, pokolenia)
     * @param racingCfg - parametry Racing
     */
    void runMetaEAWithRacing(const MetaEAConfig &metaCfg,
                             const Racing::Configuration &racingCfg)
    {
        MetaEA meta(metaCfg, racingCfg);

        // Funkcja-lambda wywołująca Twój AE z danym zestawem parametrów.
        auto evaluateFunc = [this](const ParameterSet &ps) {
            return this->runOneEvaluation(ps);
        };

        ParameterSet best = meta.runMetaEA(evaluateFunc);
        
        // Opcjonalnie można zapisać tu wynik.
        std::cout << "Best found param set: " << best.toString() << std::endl;
    }

private:
    /**
     * Tu jest kluczowy punkt integracji z Twoim kodem AE.
     * *Uruchamiamy Twój algorytm* z podanym zestawem parametrów.
     *  W zależności od stylu:
     *    - Możesz dynamicznie ustawić `GAConfig::getInstance()` zmienne
     *      (populationSize, mutationRate, itp.) i potem wywołać normalnie `runGeneticAlgorithm(...)`
     *    - Albo ładujesz nowy plik .cfg
     *    - Zwracasz TuningResult (fitness, czas, ewentualnie coverage).
     */
    TuningResult runOneEvaluation(const ParameterSet& parameterSet) {
        // Create a new instance for testing
        DNAInstance instance;
        instance.setK(10);
        instance.setDNA("ACGTACGTACGTACGTACGT");
        instance.setSpectrum({"ACGT", "CGTA", "GTAC", "TACG"});

        // Run genetic algorithm with the parameters
        double fitness = 0.0;
        runGeneticAlgorithm(instance, "output.txt", 100, "log.txt", false);

        // Return the result with fitness and execution time
        return TuningResult(parameterSet, fitness, 0.0);
    }

    /**
     * Zapis do pliku CSV wyników – z kolumnami jak w TuningResult::toCSV().
     */
    void saveResults(const std::vector<TuningResult> &results)
    {
        std::ofstream fout(m_outputFile, std::ios::app);
        if (!fout.is_open()) {
            std::cerr << "Failed to open results file: " << m_outputFile << "\n";
            return;
        }

        // Warto dodać nagłówek (jeśli plik jest świeży) – tutaj wstawiamy zawsze:
        fout << "Testowane parametry;Uzyskany fitness;Czas;Dodatkowe metryki\n";
        for (auto &r : results) {
            fout << r.toCSV() << "\n";
        }
        fout.close();
    }

    std::string m_outputFile;
};

#endif //PARAMETER_TUNING_MANAGER_H
