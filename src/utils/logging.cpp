#include "utils/logging.h"
#include <fstream>
#include <iostream>
#include <mutex>
#include <ctime>
#include <sstream>
#include <iomanip>

namespace {
    std::mutex logMutex;
    std::ofstream logFile;
    bool isInitialized = false;
    LogLevel currentLevel = LogLevel::INFO;
    
    std::string getCurrentTimestamp() {
        auto now = std::chrono::system_clock::now();
        auto time = std::chrono::system_clock::to_time_t(now);
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()) % 1000;
            
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S")
           << '.' << std::setfill('0') << std::setw(3) << ms.count();
        return ss.str();
    }
    
    std::string levelToString(LogLevel level) {
        switch (level) {
            case LogLevel::DEBUG: return "DEBUG";
            case LogLevel::INFO: return "INFO";
            case LogLevel::WARNING: return "WARNING";
            case LogLevel::ERROR: return "ERROR";
            default: return "UNKNOWN";
        }
    }
}

void Logger::initialize(const std::string& filename) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (!isInitialized) {
        logFile.open(filename, std::ios::out | std::ios::app);
        if (!logFile.is_open()) {
            std::cerr << "Failed to open log file: " << filename << std::endl;
            return;
        }
        isInitialized = true;
    }
}

void Logger::cleanup() {
    std::lock_guard<std::mutex> lock(logMutex);
    if (isInitialized && logFile.is_open()) {
        logFile.close();
        isInitialized = false;
    }
}

void Logger::setLogLevel(LogLevel level) {
    std::lock_guard<std::mutex> lock(logMutex);
    currentLevel = level;
}

void Logger::log(LogLevel level, const std::string& message, const char* file, int line) {
    if (level < currentLevel) return;
    
    std::lock_guard<std::mutex> lock(logMutex);
    if (!isInitialized) {
        std::cerr << "Logger not initialized" << std::endl;
        return;
    }
    
    std::ostringstream ss;
    ss << getCurrentTimestamp() << " ["
       << std::setw(7) << std::left << levelToString(level) << "] "
       << file << ":" << line << " - "
       << message << std::endl;
    
    std::string output = ss.str();
    logFile << output;
    logFile.flush();
    
    // Print to console for all levels in debug mode or ERROR level
    if (level == LogLevel::ERROR || currentLevel == LogLevel::DEBUG) {
        if (level == LogLevel::ERROR) {
            std::cerr << output;
        } else {
            std::cout << output;
        }
    }
}