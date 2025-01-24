#pragma once
#include <string>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <iostream>

class Logger {
private:
    static std::ofstream logFile;
    static std::streambuf* coutBuffer;
    static std::streambuf* cerrBuffer;

public:
    static void init() {
        logFile.open("sbh_framework.log", std::ios::out | std::ios::trunc);
        if (logFile.is_open()) {
            coutBuffer = std::cout.rdbuf();
            cerrBuffer = std::cerr.rdbuf();
            std::cout.rdbuf(logFile.rdbuf());
            std::cerr.rdbuf(logFile.rdbuf());
        }
    }

    static void log(const std::string& level, const std::string& message) {
        auto now = std::time(nullptr);
        auto* timeinfo = std::localtime(&now);
        char timestamp[25];
        std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo);
        logFile << timestamp << " - " << level << " - " << message << std::endl;
    }

    static void cleanup() {
        if (logFile.is_open()) {
            std::cout.rdbuf(coutBuffer);
            std::cerr.rdbuf(cerrBuffer);
            logFile.close();
        }
    }
};

// Global logging macros
#define LOG(level, msg) Logger::log(std::string(level), std::string(msg))
#define LOG_DEBUG(msg) Logger::log("DEBUG", std::string(msg))
#define LOG_INFO(msg) Logger::log("INFO", std::string(msg))
#define LOG_WARNING(msg) Logger::log("WARNING", std::string(msg))
#define LOG_ERROR(msg) Logger::log("ERROR", std::string(msg))

// For conditional debug logging
#ifdef _DEBUG
#define DEBUG_LOG(msg) LOG_DEBUG(msg)
#else
#define DEBUG_LOG(msg)
#endif 