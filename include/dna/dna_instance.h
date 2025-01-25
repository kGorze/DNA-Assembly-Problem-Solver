#pragma once

#include "utils/logging.h"
#include <string>
#include <vector>
#include <mutex>

class DNAInstance {
private:
    int n = 0;  // DNA length
    int k = 0;
    int deltaK = 0;
    int lNeg = 0;
    int lPoz = 0;
    bool repAllowed = false;
    int probablePositive = 0;
    int startIndex = -1;
    int size = 0;
    std::string targetSequence;
    std::string m_dna;
    std::string m_originalDNA;
    std::vector<std::string> m_spectrum;
    mutable std::mutex m_mutex;  // For thread-safe access to DNA and spectrum

public:
    DNAInstance() = default;
    
    DNAInstance(const std::string& originalDNA, int kValue)
        : k(std::max(1, kValue))
        , m_originalDNA(originalDNA) {}
    
    // Move constructor
    DNAInstance(DNAInstance&& other) noexcept
        : n(other.n)
        , k(other.k)
        , deltaK(other.deltaK)
        , lNeg(other.lNeg)
        , lPoz(other.lPoz)
        , repAllowed(other.repAllowed)
        , probablePositive(other.probablePositive)
        , startIndex(other.startIndex)
        , size(other.size)
        , targetSequence(std::move(other.targetSequence))
        , m_dna(std::move(other.m_dna))
        , m_originalDNA(std::move(other.m_originalDNA))
        , m_spectrum(std::move(other.m_spectrum)) {}
    
    // Move assignment operator
    DNAInstance& operator=(DNAInstance&& other) noexcept {
        if (this != &other) {
            n = other.n;
            k = other.k;
            deltaK = other.deltaK;
            lNeg = other.lNeg;
            lPoz = other.lPoz;
            repAllowed = other.repAllowed;
            probablePositive = other.probablePositive;
            startIndex = other.startIndex;
            size = other.size;
            targetSequence = std::move(other.targetSequence);
            m_dna = std::move(other.m_dna);
            m_originalDNA = std::move(other.m_originalDNA);
            m_spectrum = std::move(other.m_spectrum);
        }
        return *this;
    }
    
    // Delete copy constructor and assignment operator
    DNAInstance(const DNAInstance&) = delete;
    DNAInstance& operator=(const DNAInstance&) = delete;
    
    // Getters for instance parameters
    int getN() const { return n; }
    int getK() const { return k; }
    int getDeltaK() const { return deltaK; }
    int getLNeg() const { return lNeg; }
    int getLPoz() const { return lPoz; }
    bool isRepAllowed() const { return repAllowed; }
    int getProbablePositive() const { return probablePositive; }
    int getStartIndex() const { return startIndex; }
    int getSize() const { return size; }
    const std::string& getTargetSequence() const { return targetSequence; }
    
    // Thread-safe getters for DNA and spectrum
    std::string getDNA() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_dna; 
    }
    
    std::vector<std::string> getSpectrum() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_spectrum;
    }
    
    // Thread-safe getter for modifiable spectrum
    std::vector<std::string>& getModifiableSpectrum() {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_spectrum;
    }
    
    const std::string& getOriginalDNA() const { return m_originalDNA; }
    
    // Setters with validation
    void setN(int value) { 
        if (value <= 0) {
            LOG_WARNING("Invalid N value: " + std::to_string(value) + ", using 1");
            n = 1;
        } else {
            n = value;
        }
    }
    
    void setK(int value) {
        if (value <= 0) {
            LOG_WARNING("Invalid K value: " + std::to_string(value) + ", using 1");
            k = 1;
        } else {
            k = value;
        }
    }
    
    void setDeltaK(int value) {
        if (value < 0) {
            LOG_WARNING("Invalid deltaK value: " + std::to_string(value) + ", using 0");
            deltaK = 0;
        } else {
            deltaK = value;
        }
    }
    
    void setLNeg(int value) {
        if (value < 0) {
            LOG_WARNING("Invalid lNeg value: " + std::to_string(value) + ", using 0");
            lNeg = 0;
        } else {
            lNeg = value;
        }
    }
    
    void setLPoz(int value) {
        if (value < 0) {
            LOG_WARNING("Invalid lPoz value: " + std::to_string(value) + ", using 0");
            lPoz = 0;
        } else {
            lPoz = value;
        }
    }
    
    void setRepAllowed(bool value) { repAllowed = value; }
    
    void setProbablePositive(int value) {
        if (value < 0) {
            LOG_WARNING("Invalid probablePositive value: " + std::to_string(value) + ", using 0");
            probablePositive = 0;
        } else {
            probablePositive = value;
        }
    }
    
    void setDNA(const std::string& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (value.empty()) {
            LOG_WARNING("Empty DNA sequence provided");
        }
        m_dna = value;
    }
    
    void setSpectrum(const std::vector<std::string>& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (value.empty()) {
            LOG_WARNING("Empty spectrum provided");
        }
        m_spectrum = value;
    }
    
    void setStartIndex(int value) {
        if (value < -1) {
            LOG_WARNING("Invalid startIndex value: " + std::to_string(value) + ", using -1");
            startIndex = -1;
        } else {
            startIndex = value;
        }
    }
    
    void setSize(int value) {
        if (value <= 0) {
            LOG_WARNING("Invalid size value: " + std::to_string(value) + ", using 1");
            size = 1;
        } else {
            size = value;
        }
    }
    
    void setTargetSequence(const std::string& value) {
        if (value.empty()) {
            LOG_WARNING("Empty target sequence provided");
        }
        targetSequence = value;
    }

    // Additional functionality
    int findStartVertexIndex(const DNAInstance& instance);
}; 