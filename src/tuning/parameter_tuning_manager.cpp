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

std::vector<TuningResult> ParameterTuningManager::runRacingOnly(
    const std::vector<ParameterSet>& candidates,
    const Racing::Configuration& config,
    std::function<TuningResult(const ParameterSet&)> evaluateFunc) {
    
    // 1. Collect all keys from all candidates
    std::set<std::string> allKeys;
    for (const auto& ps : candidates) {
        for (const auto& kv : ps.params) {
            allKeys.insert(kv.first);
        }
    }
    std::vector<std::string> headerKeys(allKeys.begin(), allKeys.end());

    // 2. Open CSV file
    std::ofstream outFile(m_outputFile);
    if (!outFile.is_open()) {
        std::cerr << "Error: cannot open file " << m_outputFile << "\n";
        return {};
    }

    // 3. Write header: all parameters + fitness + executionTime
    for (const auto& key : headerKeys) {
        outFile << key << ",";
    }
    outFile << "fitness,executionTime\n";

    // 4. Run racing
    Racing::Manager racing(config);
    std::vector<TuningResult> results;
    results.reserve(candidates.size());
    
    size_t total = candidates.size();
    size_t tested = 0;

    for (const auto& candidate : candidates) {
        // Evaluate candidate
        TuningResult tr = evaluateFunc(candidate);
        results.push_back(tr);

        // Write parameter values
        for (const auto& key : headerKeys) {
            auto it = candidate.params.find(key);
            if (it != candidate.params.end()) {
                outFile << it->second;
            }
            outFile << ","; // separator
        }

        // Write fitness and executionTime
        outFile << tr.fitness << "," << tr.executionTime << "\n";
        outFile.flush();

        // Progress
        tested++;
        double percent = 100.0 * tested / total;
        std::cout << "[TUNING] Tested " << tested << " / " << total
                  << " (" << std::fixed << std::setprecision(1)
                  << percent << "%)\n";
    }

    outFile.close();
    std::cout << "Tuning finished. Results saved to: " << m_outputFile << std::endl;

    // Select best candidates
    auto selectedResults = racing.runRacing(candidates, evaluateFunc);
    return selectedResults;
}