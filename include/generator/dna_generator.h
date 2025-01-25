//
// Created by konrad_guest on 28/12/2024.
// SMART

#ifndef DNA_GENERATOR_H
#define DNA_GENERATOR_H

#include <string>
#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <mutex>
#include <stdexcept>
#include "utils/logging.h"
#include "dna/dna_instance.h"
#include "dna/dna_instance_io.h"

/**
 * Klasa generująca DNA (ciąg znaków 'ACGT') o zadanej długości n.
 * Jeżeli repAllowed = false, można tu rozszerzyć logikę ograniczającą
 * częste powtórzenia (w przykładzie pomijamy, bo zależy mocno od definicji
 * "powtórzeń").
 */
class DNAGenerator {
private:
    int n;
    int k;
    int deltaK;
    mutable std::mutex m_mutex;
    
    bool validateParameters() {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (n <= 0) {
            LOG_ERROR("Invalid n value: " + std::to_string(n));
            return false;
        }
        if (k <= 0) {
            LOG_ERROR("Invalid k value: " + std::to_string(k));
            return false;
        }
        if (deltaK < 0) {
            LOG_ERROR("Invalid deltaK value: " + std::to_string(deltaK));
            return false;
        }
        return true;
    }

public:
    DNAGenerator() : n(0), k(0), deltaK(0) {}
    
    void setParameters(int n_, int k_, int deltaK_) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (n_ <= 0) throw std::invalid_argument("n must be positive");
        if (k_ <= 0) throw std::invalid_argument("k must be positive");
        if (deltaK_ < 0) throw std::invalid_argument("deltaK cannot be negative");
        
        n = n_;
        k = k_;
        deltaK = deltaK_;
    }
    
    std::string generateDNA(int length, bool repAllowed) {
        if (length <= 0) {
            throw std::invalid_argument("DNA length must be positive");
        }
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            static const char nucleotides[] = {'A', 'C', 'G', 'T'};
            std::string dna;
            dna.reserve(length);
            
            auto& rng = RandomGenerator::getInstance().get();
            std::uniform_int_distribution<> dis(0, 3);
            
            for (int i = 0; i < length; ++i) {
                dna.push_back(nucleotides[dis(rng)]);
            }
            
            if (!repAllowed) {
                // Here you could add logic to prevent repetitions
                // For now, we'll just log that this feature is not implemented
                LOG_WARNING("Repetition prevention not implemented");
            }
            
            return dna;
        } catch (const std::exception& e) {
            LOG_ERROR("Error generating DNA: " + std::string(e.what()));
            throw;
        }
    }
    
    bool generate() {
        std::lock_guard<std::mutex> lock(m_mutex);
        return validateParameters();
    }
    
    DNAInstance generateRandomInstance(int size) {
        if (size <= 0) {
            throw std::invalid_argument("Instance size must be positive");
        }
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            DNAInstance instance;
            instance.setSize(size);
            
            std::string sequence = generateDNA(size, instance.isRepAllowed());
            instance.setDNA(sequence);
            instance.setTargetSequence(sequence);
            
            return instance;
        } catch (const std::exception& e) {
            LOG_ERROR("Error generating random instance: " + std::string(e.what()));
            throw;
        }
    }
    
    bool saveToFile(const DNAInstance& instance, const std::string& filename) {
        if (filename.empty()) {
            throw std::invalid_argument("Filename cannot be empty");
        }
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            std::ofstream file(filename);
            if (!file) {
                LOG_ERROR("Could not open file for writing: " + filename);
                return false;
            }
            
            file << instance.getSize() << "\n";
            file << instance.getTargetSequence() << "\n";
            return true;
        } catch (const std::exception& e) {
            LOG_ERROR("Error saving to file: " + std::string(e.what()));
            return false;
        }
    }
    
    DNAInstance loadFromFile(const std::string& filename) {
        if (filename.empty()) {
            throw std::invalid_argument("Filename cannot be empty");
        }
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            std::ifstream file(filename);
            if (!file) {
                throw std::runtime_error("Could not open file for reading: " + filename);
            }
            
            DNAInstance instance;
            int size;
            std::string sequence;
            
            file >> size;
            file >> sequence;
            
            if (size <= 0) {
                throw std::runtime_error("Invalid size in file: " + std::to_string(size));
            }
            
            instance.setSize(size);
            instance.setTargetSequence(sequence);
            instance.setDNA(sequence);
            
            return instance;
        } catch (const std::exception& e) {
            LOG_ERROR("Error loading from file: " + std::string(e.what()));
            throw;
        }
    }
};

/**
 * Klasa odpowiedzialna za generowanie spektrum (k-merów),
 * z uwzględnieniem reguł o zmiennej długości (deltaK).
 */
class SpectrumGenerator {
public:
    /**
     * Generuje wektor oligonukleotydów według zasad:
     *  - Pierwszy oligonukleotyd zawsze ma długość dokładnie k.
     *  - Kolejne (oprócz ostatnich k+2) mogą mieć długość k ± wylosowana_wartość (z zakresu 0..deltaK).
     *    Jeśli wylosowano > 0, to z 50% szansą plus, 50% minusem.
     *  - Ostatnich k+2 oligonukleotydów zawsze ma długość k.
     *  - Okno przesuwamy zawsze o 1 pozycję w prawo.
     */
    std::vector<std::string> generateSpectrum(const std::string &dna, int k, int deltaK);
};

/**
 * Singleton zarządzający globalnym generatorem liczb pseudolosowych.
 */
class RandomGenerator {
private:
    RandomGenerator() {
        std::random_device rd;
        gen.seed(rd() ^ (unsigned long long)std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
    
    std::mt19937 gen;
    mutable std::mutex m_mutex;

public:
    static RandomGenerator& getInstance() {
        static RandomGenerator instance;
        return instance;
    }

    std::mt19937& get() {
        std::lock_guard<std::mutex> lock(m_mutex);
        return gen;
    }

    // Delete copy constructor and assignment operator
    RandomGenerator(const RandomGenerator&) = delete;
    RandomGenerator& operator=(const RandomGenerator&) = delete;
};

/**
 * Interfejs wprowadzania błędów (negatywnych/pozytywnych).
 */
class IErrorIntroductionStrategy {
public:
    virtual ~IErrorIntroductionStrategy() = default;
    virtual void introduceErrors(DNAInstance &instance) = 0;
};

/**
 * Kontekst strategii błędów.
 */
class ErrorContext {
private:
    IErrorIntroductionStrategy* strategy = nullptr;

public:
    void setStrategy(IErrorIntroductionStrategy* s) {
        strategy = s;
    }

    void execute(DNAInstance &instance) {
        if(strategy) {
            strategy->introduceErrors(instance);
        }
    }
};

/**
 * Strategia wprowadzania błędów negatywnych:
 * - usuwa z wektora 'spectrum' lNeg losowych elementów (o ile to możliwe).
 */
class NegativeErrorIntroducer : public IErrorIntroductionStrategy {
private:
    int lNeg;
public:
    NegativeErrorIntroducer(int numNeg) : lNeg(numNeg) {}
    void introduceErrors(DNAInstance &instance) override;
};

/**
 * Strategia wprowadzania błędów pozytywnych:
 * - w zależności od probablePositive:
 *    0: dodaje lPoz losowych sekwencji o długości k ± deltaK
 *    1: duplikuje istniejące oligo w parze, modyfikując 2 nukleotydy
 */
class PositiveErrorIntroducer : public IErrorIntroductionStrategy {
private:
    int lPoz;
public:
    PositiveErrorIntroducer(int numPoz) : lPoz(numPoz) {}
    void introduceErrors(DNAInstance &instance) override;
};

#endif // DNA_GENERATOR_H
