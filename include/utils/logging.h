#pragma once

#include <string>
#include <sstream>

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

template<typename T>
void log_value(std::ostringstream& ss, const T& value) {
    ss << value;
}

template<typename... Args>
std::string format_log(const std::string& fmt, const Args&... args) {
    if constexpr (sizeof...(args) == 0) {
        return fmt;  // Return format string as-is when no arguments provided
    }
    
    std::ostringstream ss;
    size_t pos = 0;
    size_t lastPos = 0;
    
    while ((pos = fmt.find("{}", pos)) != std::string::npos) {
        ss << fmt.substr(lastPos, pos - lastPos);
        int dummy[] = {0, ((void)log_value(ss, args), 0)...};
        (void)dummy;  // Suppress unused variable warning
        lastPos = pos + 2;  // Skip "{}"
        pos += 2;
    }
    ss << fmt.substr(lastPos);
    return ss.str();
}

#define LOG_DEBUG(fmt, ...) Logger::log(LogLevel::DEBUG, format_log(fmt, ##__VA_ARGS__), __FILE__, __LINE__)
#define LOG_INFO(fmt, ...) Logger::log(LogLevel::INFO, format_log(fmt, ##__VA_ARGS__), __FILE__, __LINE__)
#define LOG_WARNING(fmt, ...) Logger::log(LogLevel::WARNING, format_log(fmt, ##__VA_ARGS__), __FILE__, __LINE__)
#define LOG_ERROR(fmt, ...) Logger::log(LogLevel::ERROR, format_log(fmt, ##__VA_ARGS__), __FILE__, __LINE__)