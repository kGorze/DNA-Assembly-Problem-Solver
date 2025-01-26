//
// Created by konrad_guest on 23/01/2025.
//
#include "../../include/tuning/parameter_tuning_manager.h"
#include "../../include/tuning/racing.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <set>
#include <algorithm>

void ParameterTuningManager::runRacingOnly(
    const std::vector<ParameterSet> &candidateParams,
    const Racing::Configuration &rc,
    std::function<TuningResult(const ParameterSet&)> evaluateFunc)
{
    // 1. Zbieramy wszystkie klucze ze wszystkich kandydatów
    std::set<std::string> allKeys;
    for (auto &ps : candidateParams) {
        for (auto &kv : ps.params) {
            allKeys.insert(kv.first);
        }
    }
    std::vector<std::string> headerKeys(allKeys.begin(), allKeys.end());

    // 2. Otwieramy plik CSV
    std::ofstream outFile(m_outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: cannot open file " << m_outputFile << "\n";
        return;
    }

    // 3. Wypisujemy nagłówek: wszystkie parametry + fitness + executionTime
    for (auto &key : headerKeys) {
        outFile << key << ",";
    }
    outFile << "fitness,executionTime\n";

    // 4. Testowanie kandydatów
    size_t total = candidateParams.size();
    size_t tested = 0;

    for (auto &ps : candidateParams) {
        // Wywołanie funkcji oceniającej
        TuningResult tr = evaluateFunc(ps);

        // Zapis wartości parametrów
        for (auto &key : headerKeys) {
            auto it = ps.params.find(key);
            if (it != ps.params.end()) {
                outFile << it->second;
            }
            outFile << ","; // separator
        }

        // Na koniec: fitness, executionTime
        outFile << tr.fitness << "," << tr.executionTime << "\n";

        outFile.flush();
        // Postęp
        tested++;
        double percent = 100.0 * tested / total;
        std::cout << "[TUNING] Tested " << tested << " / " << total
                  << " (" << std::fixed << std::setprecision(1)
                  << percent << "%)\n";
    }

    outFile.close();

    std::cout << "Tuning finished. Results saved to: " << m_outputFile << std::endl;
}

std::vector<TuningResult> ParameterTuningManager::runRacingOnly(
    const std::vector<ParameterSet>& candidates,
    const Racing::Configuration& config,
    std::function<TuningResult(const ParameterSet&)> evaluateFunc) {
    
    Racing racing(config);
    std::vector<TuningResult> results;
    results.reserve(candidates.size());
    
    for (const auto& candidate : candidates) {
        results.push_back(evaluateFunc(candidate));
    }
    
    auto selectedIndices = racing.selectBest(results);
    std::vector<TuningResult> selectedResults;
    selectedResults.reserve(selectedIndices.size());
    
    for (size_t idx : selectedIndices) {
        selectedResults.push_back(results[idx]);
    }
    
    return selectedResults;
}