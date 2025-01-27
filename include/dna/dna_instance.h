#pragma once

#include "utils/logging.h"
#include "utils/random.h"
#include <string>
#include <vector>
#include <mutex>
#include <random>
#include <memory>
#include <algorithm>

class DNAInstance {
private:
    int n = 0;  // DNA length
    int k = 0;
    int deltaK = 0;
    int lNeg = 0;
    int lPoz = 0;
    bool repAllowed = false;
    double probablePositive = 0.0;
    int startIndex = 0;
    int size = 0;
    std::string targetSequence;
    std::string m_dna;
    std::string m_originalDNA;
    std::vector<std::string> m_spectrum;
    mutable std::mutex m_mutex;  // For thread-safe access to DNA and spectrum
    std::unique_ptr<Random> m_random = std::make_unique<Random>();

    void generateSpectrum();
    std::string generateRandomDNA(int length, Random& random) const;
    bool validateSpectrum() const;

public:
    DNAInstance() = default;
    
    DNAInstance(const DNAInstance& other) {
        std::lock_guard<std::mutex> lock(other.m_mutex);
        n = other.n;
        k = other.k;
        deltaK = other.deltaK;
        lNeg = other.lNeg;
        lPoz = other.lPoz;
        repAllowed = other.repAllowed;
        probablePositive = other.probablePositive;
        startIndex = other.startIndex;
        size = other.size;
        targetSequence = other.targetSequence;
        m_dna = other.m_dna;
        m_originalDNA = other.m_originalDNA;
        m_spectrum = other.m_spectrum;
        m_random = std::make_unique<Random>();
    }
    
    // Move constructor
    DNAInstance(DNAInstance&& other) noexcept {
        std::lock_guard<std::mutex> lock(other.m_mutex);
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
        m_random = std::move(other.m_random);
    }
    
    // Copy assignment
    DNAInstance& operator=(const DNAInstance& other) {
        if (this != &other) {
            std::unique_lock<std::mutex> lock1(m_mutex, std::defer_lock);
            std::unique_lock<std::mutex> lock2(other.m_mutex, std::defer_lock);
            std::lock(lock1, lock2);  // Prevent deadlock
            
            n = other.n;
            k = other.k;
            deltaK = other.deltaK;
            lNeg = other.lNeg;
            lPoz = other.lPoz;
            repAllowed = other.repAllowed;
            probablePositive = other.probablePositive;
            startIndex = other.startIndex;
            size = other.size;
            targetSequence = other.targetSequence;
            m_dna = other.m_dna;
            m_originalDNA = other.m_originalDNA;
            m_spectrum = other.m_spectrum;
            m_random = std::make_unique<Random>();
        }
        return *this;
    }
    
    // Move assignment
    DNAInstance& operator=(DNAInstance&& other) noexcept {
        if (this != &other) {
            std::unique_lock<std::mutex> lock1(m_mutex, std::defer_lock);
            std::unique_lock<std::mutex> lock2(other.m_mutex, std::defer_lock);
            std::lock(lock1, lock2);  // Prevent deadlock
            
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
            m_random = std::move(other.m_random);
        }
        return *this;
    }

    // Getters and setters with thread safety
    int getN() const { return n; }
    void setN(int value) { n = value; }
    
    int getK() const { return k; }
    void setK(int value) { k = value; }
    
    int getDeltaK() const { return deltaK; }
    void setDeltaK(int value) { deltaK = value; }
    
    int getLNeg() const { return lNeg; }
    void setLNeg(int value) { lNeg = value; }
    
    int getLPoz() const { return lPoz; }
    void setLPoz(int value) { lPoz = value; }
    
    bool isRepAllowed() const { return repAllowed; }
    void setRepAllowed(bool value) { repAllowed = value; }
    
    double getProbablePositive() const { return probablePositive; }
    void setProbablePositive(double value) { probablePositive = value; }
    
    int getStartIndex() const { return startIndex; }
    void setStartIndex(int value) { startIndex = value; }
    
    int getSize() const { return size; }
    void setSize(int value) { size = value; }
    
    const std::string& getTargetSequence() const { return targetSequence; }
    void setTargetSequence(const std::string& value) { targetSequence = value; }
    
    const std::string& getDNA() const { return m_dna; }
    void setDNA(const std::string& value) { m_dna = value; }
    
    const std::string& getOriginalDNA() const { return m_originalDNA; }
    void setOriginalDNA(const std::string& value) { m_originalDNA = value; }
    
    const std::vector<std::string>& getSpectrum() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_spectrum;
    }
    
    void setSpectrum(const std::vector<std::string>& value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_spectrum = value;
    }
    
    void clearSpectrum() {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_spectrum.clear();
    }

    // Additional functionality
    int findStartVertexIndex(const DNAInstance& instance);

    // Constructor for test cases
    DNAInstance(int n, int k, int lNeg, int lPoz, int maxErrors, bool allowNegative, double errorProb, int seed);

    double calculateFitness(const std::string& solution) const {
        // Calculate fitness based on how well the solution matches the spectrum
        // Higher fitness means better match
        double fitness = 0.0;
        
        // Convert solution to k-mers
        std::vector<std::string> solutionKmers;
        for (size_t i = 0; i <= solution.length() - k; ++i) {
            solutionKmers.push_back(solution.substr(i, k));
        }
        
        // Count matching k-mers
        std::lock_guard<std::mutex> lock(m_mutex);
        for (const auto& kmer : solutionKmers) {
            if (std::find(m_spectrum.begin(), m_spectrum.end(), kmer) != m_spectrum.end()) {
                fitness += 1.0;
            }
        }
        
        return fitness / m_spectrum.size();  // Normalize by spectrum size
    }
}; 