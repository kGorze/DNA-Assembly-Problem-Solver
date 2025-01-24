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


/**
 * Główna klasa spinająca cały proces tuningu.
 * Można jej użyć do:
 *   - Samodzielnego trybu "tylko Racing" (podajemy listę kandydatów).
 *   - Trybu "Meta-EA + Racing" (szukanie optymalnych parametrów w dużej przestrzeni).
 */
class ParameterTuningManager {
public:
    ParameterTuningManager(const std::string &outputCsvFile)
    : m_outputFile(outputCsvFile)
    {
    }

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

    void runRacingOnly(const std::vector<ParameterSet> &candidateParameters,
                   const Racing::Configuration &racingCfg,
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
    TuningResult runOneEvaluation(const ParameterSet& params) {
        TuningResult result;
        try {
            // Create a test instance using DNAInstanceBuilder instead
            DNAInstanceBuilder builder;
            builder.setN(100)
                   .setK(7)
                   .setDeltaK(1)
                   .setLNeg(0)
                   .setLPoz(0)
                   .setRepAllowed(true)
                   .buildDNA()
                   .buildSpectrum();
            
            DNAInstance instance = builder.getInstance();
            
            // Run GA with these parameters
            auto start = std::chrono::high_resolution_clock::now();
            double fitness = runGeneticAlgorithmWrapper(instance);
            auto end = std::chrono::high_resolution_clock::now();
            
            result.fitness = fitness;
            result.executionTime = std::chrono::duration<double>(end - start).count();
            result.parameterSet = params;
            
        } catch (const std::exception& e) {
            result.fitness = 0.0;  // or some other invalid fitness value
            result.executionTime = 0.0;
            result.parameterSet = params;
        }
        return result;
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
