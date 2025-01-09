//
// Created by konrad_guest on 09/01/2025.
//

#include "utils/performance_profilling_framework.h"


// Initialize static member
Profiler* Profiler::instance = nullptr;

Profiler::Profiler() : programStart(std::chrono::steady_clock::now()) {}

Profiler& Profiler::getInstance() {
    if (!instance) {
        instance = new Profiler();
    }
    return *instance;
}

Profiler::ScopedTimer::ScopedTimer(const std::string& name)
    : functionName(name)
    , start(std::chrono::steady_clock::now())
    , profiler(Profiler::getInstance()) {}

Profiler::ScopedTimer::~ScopedTimer() {
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    profiler.addMeasurement(functionName, duration);
}

void Profiler::addMeasurement(const std::string& functionName, std::chrono::nanoseconds duration) {
    std::lock_guard<std::mutex> lock(statsMutex);
    auto& funcStats = stats[functionName];
    funcStats.functionName = functionName;
    funcStats.totalTime += duration;
    funcStats.callCount++;
}

void Profiler::saveReport(const std::string& filename) {
    auto programEnd = std::chrono::steady_clock::now();
    auto totalProgramTime = std::chrono::duration_cast<std::chrono::nanoseconds>(
        programEnd - programStart).count();

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // CSV header
    file << "Function Name,Total Time (ms),Calls,Avg Time per Call (ms),Percentage of Total Time\n";

    // Convert stats to vector for sorting
    std::vector<FunctionStats> sortedStats;
    for (const auto& pair : stats) {
        sortedStats.push_back(pair.second);
    }

    // Sort by total time (descending)
    std::sort(sortedStats.begin(), sortedStats.end(),
        [](const FunctionStats& a, const FunctionStats& b) {
            return a.totalTime > b.totalTime;
        });

    // Write stats
    for (const auto& stat : sortedStats) {
        double totalMs = stat.totalTime.count() / 1e6; // Convert ns to ms
        double avgMs = totalMs / stat.callCount;
        double percentage = (stat.totalTime.count() * 100.0) / totalProgramTime;
        
        file << stat.functionName << ","
             << totalMs << ","
             << stat.callCount << ","
             << avgMs << ","
             << percentage << "\n";
    }

    file.close();
    std::cout << "Profiling results saved to: " << filename << std::endl;
}

Profiler::~Profiler() {
    saveReport();
}