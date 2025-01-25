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
#include "dna/error_introduction.h"

/**
 * Klasa generująca DNA (ciąg znaków 'ACGT') o zadanej długości n.
 * Jeżeli repAllowed = false, można tu rozszerzyć logikę ograniczającą
 * częste powtórzenia (w przykładzie pomijamy, bo zależy mocno od definicji
 * "powtórzeń").
 */
class DNAGenerator {
private:
    mutable std::mutex m_mutex;
    int m_n = 0;
    int m_k = 0;
    int m_deltaK = 0;
    
    bool validateParameters() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (m_n <= 0) {
            LOG_ERROR("Invalid n value: " + std::to_string(m_n));
            return false;
        }
        if (m_k <= 0) {
            LOG_ERROR("Invalid k value: " + std::to_string(m_k));
            return false;
        }
        if (m_deltaK < 0) {
            LOG_ERROR("Invalid deltaK value: " + std::to_string(m_deltaK));
            return false;
        }
        return true;
    }

public:
    DNAGenerator() = default;
    
    void setParameters(int n, int k, int deltaK) {
        std::lock_guard<std::mutex> lock(m_mutex);
        if (n <= 0) throw std::invalid_argument("n must be positive");
        if (k <= 0) throw std::invalid_argument("k must be positive");
        if (deltaK < 0) throw std::invalid_argument("deltaK cannot be negative");
        
        m_n = n;
        m_k = k;
        m_deltaK = deltaK;
    }
    
    std::string generateDNA(int length, bool repAllowed = true) const {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            if (length <= 0) {
                throw std::invalid_argument("DNA length must be positive");
            }
            
            static thread_local std::random_device rd;
            static thread_local std::mt19937 gen(rd());
            std::uniform_int_distribution<> dis(0, 3);
            
            const char nucleotides[] = {'A', 'C', 'G', 'T'};
            std::string dna;
            dna.reserve(length);
            
            for (int i = 0; i < length; ++i) {
                dna.push_back(nucleotides[dis(gen)]);
            }
            
            if (!repAllowed) {
                LOG_WARNING("Repetition prevention not implemented");
            }
            
            return dna;
        } catch (const std::exception& e) {
            LOG_ERROR("Error generating DNA: " + std::string(e.what()));
            throw;
        }
    }
    
    DNAInstance generateRandomInstance(int size, int lNeg = 0, int lPoz = 0) const {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            if (size <= 0) {
                throw std::invalid_argument("Instance size must be positive");
            }
            if (!validateParameters()) {
                throw std::runtime_error("Invalid generator parameters");
            }
            
            DNAInstance instance;
            instance.setN(size);
            instance.setK(m_k);
            instance.setDeltaK(m_deltaK);
            instance.setLNeg(lNeg);
            instance.setLPoz(lPoz);
            
            // Generate DNA and spectrum
            std::string dna = generateDNA(size);
            instance.setDNA(dna);
            
            SpectrumGenerator specGen;
            auto spectrum = specGen.generateSpectrum(dna, m_k, m_deltaK);
            instance.setSpectrum(spectrum);
            
            // Find start index before introducing errors
            std::string startFrag = dna.substr(0, m_k);
            int startIdx = -1;
            for (int i = 0; i < static_cast<int>(spectrum.size()); i++) {
                if (spectrum[i] == startFrag) {
                    startIdx = i;
                    break;
                }
            }
            instance.setStartIndex(startIdx);
            
            // Introduce errors if requested
            if (lNeg > 0) {
                auto negErr = ErrorIntroducerFactory::createNegativeErrorIntroducer(lNeg);
                negErr->introduceErrors(instance);
            }
            
            if (lPoz > 0) {
                auto posErr = ErrorIntroducerFactory::createPositiveErrorIntroducer(lPoz);
                posErr->introduceErrors(instance);
            }
            
            return instance;
        } catch (const std::exception& e) {
            LOG_ERROR("Error generating random instance: " + std::string(e.what()));
            throw;
        }
    }
    
    bool saveToFile(const DNAInstance& instance, const std::string& filename) const {
        if (filename.empty()) {
            throw std::invalid_argument("Filename cannot be empty");
        }
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            return InstanceIO::saveInstance(instance, filename);
        } catch (const std::exception& e) {
            LOG_ERROR("Error saving instance to file: " + std::string(e.what()));
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
