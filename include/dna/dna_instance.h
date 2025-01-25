#pragma once

#include "utils/logging.h"
#include <string>
#include <vector>
#include <mutex>

class DNAInstance {
public:
    DNAInstance() = default;
    
    DNAInstance(const std::string& originalDNA, int kValue)
        : m_originalDNA(originalDNA), k(std::max(1, kValue)) {}
    
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
        return dna; 
    }
    
    std::vector<std::string> getSpectrum() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return spectrum;
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
        dna = value;
    }
    
    void setSpectrum(const std::vector<std::string>& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (value.empty()) {
            LOG_WARNING("Empty spectrum provided");
        }
        spectrum = value;
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
    
    std::string dna;
    std::string targetSequence;
    std::vector<std::string> spectrum;
    std::string m_originalDNA;
    
    mutable std::mutex m_mutex;  // For thread-safe access to DNA and spectrum
}; 