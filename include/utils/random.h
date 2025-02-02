//test already in the suite.
#pragma once

#include <random>
#include <chrono>
#include <mutex>
#include <stdexcept>

class Random {
private:
    std::mt19937 m_generator;
    static constexpr size_t MAX_SIZE_T = 1000000;  // Reasonable limit for our use case

public:
    Random() {
        std::random_device rd;
        m_generator.seed(rd() ^ static_cast<unsigned long long>(std::chrono::high_resolution_clock::now().time_since_epoch().count()));
    }

    explicit Random(int seed) : m_generator(seed) {}
    
    virtual ~Random() = default;

    static Random& instance() {
        static Random instance;
        return instance;
    }

    std::mt19937& getGenerator() {
        return m_generator;
    }

    virtual double generateProbability() {
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        return dist(m_generator);
    }

    virtual int getRandomInt(int min, int max) {
        if (min > max) {
            throw std::invalid_argument("min cannot be greater than max");
        }
        std::uniform_int_distribution<int> dist(min, max);
        return dist(m_generator);
    }

    virtual size_t getRandomSizeT(size_t min, size_t max) {
        if (min > max) {
            throw std::invalid_argument("min cannot be greater than max");
        }
        if (max > MAX_SIZE_T) {
            throw std::invalid_argument("max value exceeds reasonable limit");
        }
        std::uniform_int_distribution<size_t> dist(min, max);
        return dist(m_generator);
    }
}; 