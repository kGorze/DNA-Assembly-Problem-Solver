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
    static bool initialized;

public:
    static void init() {
        if (initialized) return;
        
        logFile.open("sbh_framework.log", std::ios::out | std::ios::trunc);
        if (logFile.is_open()) {
            coutBuffer = std::cout.rdbuf();
            cerrBuffer = std::cerr.rdbuf();
            std::cout.rdbuf(logFile.rdbuf());
            std::cerr.rdbuf(logFile.rdbuf());
            initialized = true;
        }
    }

    static void log(const std::string& level, const std::string& message) {
        if (!initialized || !logFile.is_open()) return;

        try {
            auto now = std::time(nullptr);
            if (auto* timeinfo = std::localtime(&now)) {
                char timestamp[26] = {0}; // One extra byte for null terminator
                if (std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo) > 0) {
                    logFile << timestamp << " - " << level << " - " << message << std::endl;
                    logFile.flush();
                }
            }
        } catch (const std::exception&) {
            // Silently fail if logging fails - we don't want logging errors to crash the program
        }
    }

    static void cleanup() {
        if (!initialized) return;
        
        if (logFile.is_open()) {
            std::cout.rdbuf(coutBuffer);
            std::cerr.rdbuf(cerrBuffer);
            logFile.close();
        }
        initialized = false;
        coutBuffer = nullptr;
        cerrBuffer = nullptr;
    }
};

// Initialize static members
std::ofstream Logger::logFile;
std::streambuf* Logger::coutBuffer = nullptr;
std::streambuf* Logger::cerrBuffer = nullptr;
bool Logger::initialized = false;

// Global logging macros - avoid unnecessary string construction
#define LOG(level, msg) Logger::log(level, msg)
#define LOG_DEBUG(msg) Logger::log("DEBUG", msg)
#define LOG_INFO(msg) Logger::log("INFO", msg)
#define LOG_WARNING(msg) Logger::log("WARNING", msg)
#define LOG_ERROR(msg) Logger::log("ERROR", msg)

// For conditional debug logging
#ifdef _DEBUG
#define DEBUG_LOG(msg) LOG_DEBUG(msg)
#else
#define DEBUG_LOG(msg)
#endif 