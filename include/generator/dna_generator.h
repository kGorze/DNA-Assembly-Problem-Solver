//
// Created by konrad_guest on 28/12/2024.
// SMART

#pragma once

#include <string>
#include <vector>
#include <random>
#include <mutex>
#include <memory>
#include "utils/logging.h"
#include "dna/dna_instance.h"
#include "dna/dna_instance_io.h"
#include "dna/error_introduction.h"
#include "utils/random.h"

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
    
    bool validateParameters() const;

public:
    explicit DNAGenerator(std::unique_ptr<Random> random = std::make_unique<Random>());
    
    void setParameters(int n, int k, int deltaK);
    std::string generateDNA(int length, bool repAllowed = true) const;
    DNAInstance generateRandomInstance(
        int size,
        int k,
        int lNeg = 0,
        int lPoz = 0,
        int maxErrors = 0,
        bool allowNegative = false,
        double errorProb = 0.0
    ) const;
    bool saveToFile(const DNAInstance& instance, const std::string& filename) const;
    static DNAInstance loadFromFile(const std::string& filename);
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
    std::vector<std::string> generateSpectrum(const std::string& dna, int k, int deltaK);
};

/**
 * Singleton zarządzający globalnym generatorem liczb pseudolosowych.
 */
class RandomGenerator {
private:
    RandomGenerator();
    std::mt19937 gen;
    mutable std::mutex m_mutex;

public:
    static RandomGenerator& getInstance();
    std::mt19937& get();

    RandomGenerator(const RandomGenerator&) = delete;
    RandomGenerator& operator=(const RandomGenerator&) = delete;
};
