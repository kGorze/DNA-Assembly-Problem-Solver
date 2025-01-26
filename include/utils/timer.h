#pragma once

#include <chrono>
#include <string>

class Timer {
public:
    Timer() : m_start(std::chrono::high_resolution_clock::now()) {}

    void reset() {
        m_start = std::chrono::high_resolution_clock::now();
    }

    double getElapsedSeconds() const {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - m_start);
        return duration.count() / 1000000.0;
    }

    double getElapsedMilliseconds() const {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - m_start);
        return duration.count() / 1000.0;
    }

    double getElapsedMicroseconds() const {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now - m_start);
        return static_cast<double>(duration.count());
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
}; 