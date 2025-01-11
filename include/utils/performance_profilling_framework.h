//
// Created by konrad_guest on 09/01/2025.
// SMART

#ifndef PERFORMANCE_PROFILLING_FRAMEWORK_H
#define PERFORMANCE_PROFILLING_FRAMEWORK_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <mutex>
#include <iostream>

class Profiler {
private:
    struct FunctionStats {
        std::string functionName;
        std::chrono::nanoseconds totalTime{0};
        int callCount{0};
    };
    
    std::unordered_map<std::string, FunctionStats> stats;
    std::chrono::steady_clock::time_point programStart;
    std::mutex statsMutex;
    static Profiler* instance;
    
    Profiler();

public:
    static Profiler& getInstance();

    class ScopedTimer {
    private:
        std::string functionName;
        std::chrono::steady_clock::time_point start;
        Profiler& profiler;
        
    public:
        explicit ScopedTimer(const std::string& name);
        ~ScopedTimer();
    };

    void addMeasurement(const std::string& functionName, std::chrono::nanoseconds duration);
    void saveReport(const std::string& filename = "profiling_results.csv");
    ~Profiler();
};

#define PROFILE_FUNCTION() Profiler::ScopedTimer timer(__FUNCTION__)
#define PROFILE_SCOPE(name) Profiler::ScopedTimer timer(name)

#endif //PERFORMANCE_PROFILLING_FRAMEWORK_H
