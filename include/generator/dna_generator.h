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
    
    bool validateParameters() {
        return (n > 0 && k > 0 && deltaK >= 0);
    }

public:
    DNAGenerator() : n(0), k(0), deltaK(0) {}
    void setParameters(int n_, int k_, int deltaK_) {
        n = n_;
        k = k_;
        deltaK = deltaK_;
    }
    
    std::string generateDNA(int n, bool repAllowed);
    bool generate();
    DNAInstance generateRandomInstance(int size);
    bool saveToFile(const DNAInstance& instance, const std::string& filename);
    DNAInstance loadFromFile(const std::string& filename);
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
    std::mt19937 gen;

    RandomGenerator();
public:
    static RandomGenerator& getInstance();
    std::mt19937& get();

    // usuń kopiowanie
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
