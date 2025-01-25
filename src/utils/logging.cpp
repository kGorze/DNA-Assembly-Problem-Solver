#include "../../include/utils/logging.h"
#include <ctime>
#include <iomanip>
#include <filesystem>

// Define static members
std::ofstream Logger::logFile;
std::streambuf* Logger::coutBuffer = nullptr;
std::streambuf* Logger::cerrBuffer = nullptr;
bool Logger::initialized = false;
std::mutex Logger::logMutex;

void Logger::init(const std::string& logFilePath) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (initialized) return;
    
    logFile.open(logFilePath, std::ios::out | std::ios::trunc);
    if (logFile.is_open()) {
        coutBuffer = std::cout.rdbuf();
        cerrBuffer = std::cerr.rdbuf();
        std::cout.rdbuf(logFile.rdbuf());
        std::cerr.rdbuf(logFile.rdbuf());
        initialized = true;
    }
}

void Logger::cleanup() {
    std::lock_guard<std::mutex> lock(logMutex);
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

void Logger::log(const std::string& level, const std::string& message) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (!initialized || !logFile.is_open()) return;

    try {
        auto now = std::time(nullptr);
        if (auto* timeinfo = std::localtime(&now)) {
            char timestamp[26] = {0};
            if (std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo) > 0) {
                logFile << timestamp << " - " << level << " - " << message << std::endl;
                logFile.flush();
            }
        }
    } catch (const std::exception&) {
        // Silently fail if logging fails
    }
}

void Logger::log(const std::string& message, const std::string& category, 
                const std::string& file, int line, const std::string& function) {
    std::lock_guard<std::mutex> lock(logMutex);
    if (!initialized || !logFile.is_open()) return;

    try {
        auto now = std::time(nullptr);
        if (auto* timeinfo = std::localtime(&now)) {
            char timestamp[26] = {0};
            if (std::strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", timeinfo) > 0) {
                logFile << timestamp << " - " << category << " - " 
                       << file << ":" << line << " - " 
                       << function << " - " << message << std::endl;
                logFile.flush();
            }
        }
    } catch (const std::exception&) {
        // Silently fail if logging fails
    }
}

void Logger::logDebug(const std::string& message) {
    log("DEBUG", message);
}

void Logger::logInfo(const std::string& message) {
    log("INFO", message);
}

void Logger::logWarning(const std::string& message) {
    log("WARNING", message);
}

void Logger::logError(const std::string& message) {
    log("ERROR", message);
}