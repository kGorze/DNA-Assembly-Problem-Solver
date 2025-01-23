//
// Created by konrad_guest on 23/01/2025.
//

#ifndef RACING_H
#define RACING_H

#include "tuning_structures.h"
#include "../generator/dna_generator.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <functional>
#include <map>
#include <cmath>
#include <numeric>  // [ADDED LINES FOR EXTENDED RACING]

namespace Racing {

struct Configuration {
    // Parametry "wyścigu":
    double significanceLevel = 0.05;   
    int maxTrialsPerCandidate = 50;    
    int minTrialsBeforeElimination = 5; 

    // [ADDED LINES FOR EXTENDED RACING]
    bool useBootstrap = false;       // czy używamy testu bootstrap
    int bootstrapSamples = 2000;     // ile replik w bootstrapie
};

/**
 * Klasa implementująca logikę Racing do szybkiego porównywania wielu kandydatów.
 * Wywołuje "runEA(parameters)" w pętli, zbiera wyniki i eliminuje słabszych.
 */
class Manager {
public:
    Manager(const Configuration &cfg) : config(cfg) {}

    /**
     * @param candidates - lista możliwych konfiguracji parametrów
     * @param evaluateFunc - funkcja, która uruchamia Twój algorytm z danym zestawem parametrów
     *                       i zwraca osiągnięty fitness (lub inny kryterium).
     *                       Powinna także ewentualnie zwracać czas wykonania.
     */
    std::vector<TuningResult> runRacing(
            std::vector<ParameterSet> candidates,
            std::function<TuningResult(const ParameterSet&)> evaluateFunc)
    {
        // Implementacja RacingManager::runRacing
        std::vector<bool> active(candidates.size(), true);
        
        std::vector<std::vector<double>> samples(candidates.size());
        std::vector<double> bestSoFar(candidates.size(), -1e9);
        std::vector<double> totalTime(candidates.size(), 0.0);

        std::mt19937 rng(std::random_device{}());

        int aliveCount = static_cast<int>(candidates.size());
        int trial = 0;
        
        while (aliveCount > 1 && trial < config.maxTrialsPerCandidate) {
            for (size_t i = 0; i < candidates.size(); i++) {
                if (!active[i]) continue;

                TuningResult tr = evaluateFunc(candidates[i]);
                
                double fit = tr.fitness;
                double time = tr.executionTime;
                
                samples[i].push_back(fit);
                if (fit > bestSoFar[i]) {
                    bestSoFar[i] = fit;
                }
                totalTime[i] += time;
            }

            if (trial >= config.minTrialsBeforeElimination) {
                size_t bestIdx = 0;
                double bestAvg = -1e9;
                for (size_t i = 0; i < candidates.size(); i++) {
                    if (!active[i]) continue;
                    double avgFit = average(samples[i]);
                    if (avgFit > bestAvg) {
                        bestAvg = avgFit;
                        bestIdx = i;
                    }
                }

                double bestCandidateAvg = average(samples[bestIdx]);

                for (size_t i = 0; i < candidates.size(); i++) {
                    if (i == bestIdx || !active[i]) continue;
                    
                    double avgFitI = average(samples[i]);
                    
                    // [ADDED LINES FOR EXTENDED RACING]
                    // Wybieramy metodę testu: bootstrap lub Welch
                    double pValue = 1.0;
                    if (config.useBootstrap) {
                        pValue = computeBootstrapPValue(samples[i], samples[bestIdx], config.bootstrapSamples);
                    } else {
                        pValue = computeWelchTTestPValue(samples[i], samples[bestIdx]);
                    }

                    bool isMuchWorse = false;
                    if (pValue < config.significanceLevel && avgFitI < bestCandidateAvg) {
                        isMuchWorse = true;
                    }
                    
                    // Stary próg "avgFitI + 2.0 < bestCandidateAvg" (pozostawiony)
                    if (avgFitI + 2.0 < bestCandidateAvg) {
                        isMuchWorse = true;
                    }

                    if (isMuchWorse) {
                        active[i] = false;
                        aliveCount--;
                    }
                }
            }

            trial++;
        }

        std::vector<TuningResult> finalResults;
        finalResults.reserve(candidates.size());
        for (size_t i = 0; i < candidates.size(); i++) {
            if (!active[i]) continue;

            TuningResult tres;
            tres.parameterSet = candidates[i];
            tres.fitness = average(samples[i]);       
            tres.executionTime = totalTime[i];        
            finalResults.push_back(tres);
        }

        return finalResults;
    }

private:
    Configuration config;

    // Pomocnicze funkcje
    static double average(const std::vector<double> &v) {
        if (v.empty()) return 0.0;
        double sum = 0.0;
        for (auto val : v) sum += val;
        return sum / v.size();
    }

    static double stdev(const std::vector<double> &v) {
        double avg = average(v);
        double sum = 0.0;
        for (auto val : v) {
            sum += (val - avg)*(val - avg);
        }
        return (v.size() > 1) ? std::sqrt(sum / (v.size()-1)) : 0.0;
    }

    // ---------------------- WELCH T-TEST (zostaje) ----------------------
    static double computeWelchTTestPValue(const std::vector<double> &s1,
                                          const std::vector<double> &s2)
    {
        if (s1.size() < 2 || s2.size() < 2) {
            return 1.0;
        }

        double mean1 = average(s1);
        double mean2 = average(s2);
        double var1  = variance(s1);
        double var2  = variance(s2);

        double n1 = static_cast<double>(s1.size());
        double n2 = static_cast<double>(s2.size());

        double numerator = mean1 - mean2;
        double denom = std::sqrt(var1/n1 + var2/n2);
        if (denom < 1e-12) {
            return 1.0;
        }
        double tValue = numerator / denom;

        double dfNumerator = std::pow(var1/n1 + var2/n2, 2);
        double dfDenominator = std::pow(var1/n1, 2)/(n1-1.0)
                               + std::pow(var2/n2, 2)/(n2-1.0);
        double df = dfNumerator / dfDenominator;
        if (df < 1.0) df = 1.0;

        // W pełnej wersji należałoby policzyć p-value z rozkładu t.
        // Tutaj jednak w oryginale była uproszczona linia. Zostawiamy obecną, ale dopisujemy komentarz:
        // "Należy w realnym kodzie skorzystać z biblioteki stat. do cdf t-Studenta."
        double approxP = std::exp(-std::fabs(tValue));

        return approxP;
    }

    // ---------------------- BOOTSTRAP TEST ----------------------
    // Dwustronny test na różnicę średnich
    static double computeBootstrapPValue(const std::vector<double> &s1,
                                         const std::vector<double> &s2,
                                         int B)
    {
        if (s1.size() < 2 || s2.size() < 2) {
            return 1.0;
        }
        std::mt19937 rng(std::random_device{}());

        double mean1 = average(s1);
        double mean2 = average(s2);
        double observedDiff = mean1 - mean2;

        // Połączamy próbki
        std::vector<double> combined;
        combined.reserve(s1.size() + s2.size());
        combined.insert(combined.end(), s1.begin(), s1.end());
        combined.insert(combined.end(), s2.begin(), s2.end());

        double countMoreExtreme = 0.0;

        // Generujemy B replik
        for (int b = 0; b < B; b++) {
            // losujemy z U(0,1)
            std::vector<double> bootS1, bootS2;
            bootS1.reserve(s1.size());
            bootS2.reserve(s2.size());

            // "resampling with replacement" – ale z podziałem na s1.size() i s2.size()
            std::uniform_int_distribution<int> dist(0, combined.size()-1);
            for (size_t i = 0; i < s1.size(); i++) {
                int idx = dist(rng);
                bootS1.push_back(combined[idx]);
            }
            for (size_t i = 0; i < s2.size(); i++) {
                int idx = dist(rng);
                bootS2.push_back(combined[idx]);
            }

            double m1 = average(bootS1);
            double m2 = average(bootS2);
            double diff = m1 - m2;

            if (std::fabs(diff) >= std::fabs(observedDiff)) {
                countMoreExtreme += 1.0;
            }
        }

        double pValue = (countMoreExtreme + 1.0) / (B + 1.0);
        return pValue;
    }

    static double variance(const std::vector<double> &v) {
        if (v.size() < 2) return 0.0;
        double avg = average(v);
        double sumSq = 0.0;
        for (auto val : v) {
            double diff = (val - avg);
            sumSq += diff*diff;
        }
        return sumSq / (v.size() - 1);
    }
};

}; // namespace Racing

#endif //RACING_H
