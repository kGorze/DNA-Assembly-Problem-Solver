#pragma once
#include <string>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mutex>

class Logger {
public:
    static void init(const std::string& logFilePath);
    static void cleanup();
    
    static void log(const std::string& level, const std::string& message);
    static void log(const std::string& message, const std::string& category, 
                   const std::string& file, int line, const std::string& function);

    static void logDebug(const std::string& message);
    static void logInfo(const std::string& message);
    static void logWarning(const std::string& message);
    static void logError(const std::string& message);

private:
    static std::ofstream logFile;
    static std::streambuf* coutBuffer;
    static std::streambuf* cerrBuffer;
    static bool initialized;
    static std::mutex logMutex;
};

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