#pragma once

#include <string>

enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
};

class Logger {
public:
    static void initialize(const std::string& filename);
    static void cleanup();
    static void setLogLevel(LogLevel level);
    static void log(LogLevel level, const std::string& message, const char* file, int line);

private:
    Logger() = delete;
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;
};

#define LOG_DEBUG(msg) Logger::log(LogLevel::DEBUG, msg, __FILE__, __LINE__)
#define LOG_INFO(msg) Logger::log(LogLevel::INFO, msg, __FILE__, __LINE__)
#define LOG_WARNING(msg) Logger::log(LogLevel::WARNING, msg, __FILE__, __LINE__)
#define LOG_ERROR(msg) Logger::log(LogLevel::ERROR, msg, __FILE__, __LINE__) 