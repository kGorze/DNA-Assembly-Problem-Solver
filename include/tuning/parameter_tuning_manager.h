//
// Created by konrad_guest on 23/01/2025.
//

#ifndef PARAMETER_TUNING_MANAGER_H
#define PARAMETER_TUNING_MANAGER_H


#include "racing.h"
#include "parameters_parser.h"
#include "tuning_structures.h"
#include "meta_ea.h"
#include "metaheuristics/genetic_algorithm.h"
#include "../generator/dna_generator.h"


#include <fstream>
#include <functional>
#include <vector>
#include <string>
#include <chrono>



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
    TuningResult runOneEvaluation(const ParameterSet &ps)
    {
        // 1. Parametry -> GAConfig (lub inny Twój config)
        auto &cfg = GAConfig::getInstance();
        if (ps.params.find("populationSize") != ps.params.end()) {
            cfg.populationSize = std::stoi(ps.params.at("populationSize"));
        }
        if (ps.params.find("mutationRate") != ps.params.end()) {
            cfg.mutationRate = std::stod(ps.params.at("mutationRate"));
        }
        if (ps.params.find("selectionMethod") != ps.params.end()) {
            cfg.selectionMethod = ps.params.at("selectionMethod");
        }
        // Można obsłużyć tutaj inne parametry (np. crossoverRate, itp.)

        // 2. Uruchamiamy Twój algorytm:
        auto start = std::chrono::high_resolution_clock::now();

        // Tworzymy instancję DNA (lub wczytujemy z pliku), np.:
        DNAInstance instance;
        // Przykładowo wypełniamy 10 genów losowymi liczbami:
        instance.genes.resize(10, 0.0);
        {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<double> dist(0.0, 1.0);
            for (auto &g : instance.genes) {
                g = dist(rng);
            }
        }

        // Uruchamiamy GA:
        double finalFitness = runGeneticAlgorithmWrapper(instance);

        auto end = std::chrono::high_resolution_clock::now();
        double durationSec = std::chrono::duration<double>(end - start).count();

        // 3. Zwracamy wynik w strukturze TuningResult
        TuningResult tr;
        tr.parameterSet = ps;
        tr.fitness = finalFitness;
        tr.executionTime = durationSec;

        // ewentualnie wypełnić "dodatkowe metryki"
        tr.extraMetrics["Coverage"] = 50.0; // przykładowa stała
        tr.extraMetrics["EdgeScore"] = 10.0; // j.w.

        return tr;
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
