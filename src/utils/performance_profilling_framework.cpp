//
// Created by konrad_guest on 09/01/2025.
// SMART
#include "../../include/utils/performance_profilling_framework.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>

// Initialize static member
Profiler* Profiler::instance = nullptr;

Profiler::Profiler() : programStart(std::chrono::steady_clock::now()) {}

Profiler& Profiler::getInstance() {
    if (!instance) {
        instance = new Profiler();
    }
    return *instance;
}

Profiler::~Profiler() {
    try {
        saveReport();
    } catch (const std::exception& e) {
        LOG_ERROR("Error saving profiling report: " + std::string(e.what()));
    }
    delete instance;
    instance = nullptr;
}

void Profiler::addMeasurement(const std::string& functionName, std::chrono::nanoseconds duration) {
    std::lock_guard<std::mutex> lock(statsMutex);
    auto& stat = stats[functionName];
    stat.functionName = functionName;
    stat.totalTime += duration;
    stat.callCount++;
}

void Profiler::saveReport(const std::string& filename) {
    std::lock_guard<std::mutex> lock(statsMutex);
    std::ofstream file(filename);
    if (!file.is_open()) {
        LOG_ERROR("Failed to open profiling report file: " + filename);
        return;
    }

    // Calculate total program time
    auto programEnd = std::chrono::steady_clock::now();
    auto totalProgramTime = std::chrono::duration_cast<std::chrono::milliseconds>(programEnd - programStart);

    // Write header
    file << "Function Name,Total Time (ms),Call Count,Average Time (ms),% of Total Time\n";

    // Convert stats to vector for sorting
    std::vector<std::pair<std::string, FunctionStats>> sortedStats;
    for (const auto& [name, stat] : stats) {
        sortedStats.push_back({name, stat});
    }

    // Sort by total time
    std::sort(sortedStats.begin(), sortedStats.end(),
        [](const auto& a, const auto& b) {
            return a.second.totalTime > b.second.totalTime;
        });

    // Write stats
    for (const auto& [name, stat] : sortedStats) {
        auto totalTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(stat.totalTime).count();
        auto avgTimeMs = totalTimeMs / static_cast<double>(stat.callCount);
        auto percentageOfTotal = (totalTimeMs / static_cast<double>(totalProgramTime.count())) * 100.0;

        file << stat.functionName << ","
             << totalTimeMs << ","
             << stat.callCount << ","
             << std::fixed << std::setprecision(3) << avgTimeMs << ","
             << std::fixed << std::setprecision(2) << percentageOfTotal << "\n";
    }

    file.close();
    LOG_INFO("Profiling report saved to: " + filename);
}

Profiler::ScopedTimer::ScopedTimer(const std::string& name)
    : functionName(name), start(std::chrono::steady_clock::now()), profiler(Profiler::getInstance()) {}

Profiler::ScopedTimer::~ScopedTimer() {
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    profiler.addMeasurement(functionName, duration);
}
