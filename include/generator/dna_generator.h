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

/**
 * Struktura przechowująca całą instancję: DNA, parametry i wygenerowane spektrum.
 */
class DNAInstance {
private:
    std::string dna;
    int n;                // długość DNA
    int k;                // długość oligo
    int deltaK;           // maksymalna zmiana długości oligo
    int lNeg;             // liczba błędów negatywnych
    int lPoz;             // liczba błędów pozytywnych
    bool repAllowed;      // czy dozwolone powtórzenia w DNA
    int probablePositive; // sposób generowania błędów pozytywnych
    int startIndex;       // Nowy atrybut
    std::vector<std::string> spectrum;
    std::string targetSequence;
    int size;

public:
    DNAInstance() : n(0), k(0), deltaK(0), lNeg(0), lPoz(0), 
                   repAllowed(false), probablePositive(0), 
                   startIndex(-1), size(0) {}

    // Add these getters/setters
    const std::string& getTargetSequence() const { return targetSequence; }
    void setTargetSequence(const std::string& seq) { targetSequence = seq; }
    void setSize(int s) { size = s; }
    int getSize() const { return size; }

    // get/set

    int getStartIndex() const { return startIndex; }
    void setStartIndex(int index) { startIndex = index; }
    
    void setDNA(const std::string &d) { dna = d; }
    std::string getDNA() const { return dna; }

    void setN(int val) { n = val; }
    int getN() const { return n; }

    void setK(int val) {
        if (val <= 0) {
            throw std::invalid_argument("K must be positive");
        }
        k = val;
        
        // Validate existing spectrum if any
        for (const auto& fragment : spectrum) {
            if (fragment.size() < k - deltaK || fragment.size() > k + deltaK) {
                std::string msg = "Fragment size " + std::to_string(fragment.size()) + 
                                " outside allowed range [" + std::to_string(k-deltaK) + 
                                "," + std::to_string(k+deltaK) + "]";
                LOG_WARNING(msg);
            }
        }
    }
    int getK() const { return k; }

    void setDeltaK(int val) { deltaK = val; }
    int getDeltaK() const { return deltaK; }

    void setLNeg(int val) { lNeg = val; }
    int getLNeg() const { return lNeg; }

    void setLPoz(int val) { lPoz = val; }
    int getLPoz() const { return lPoz; }

    void setRepAllowed(bool val) { repAllowed = val; }
    bool isRepAllowed() const { return repAllowed; }

    void setProbablePositive(int val) { probablePositive = val; }
    int getProbablePositive() const { return probablePositive; }
    

    void setSpectrum(const std::vector<std::string> &s) { spectrum = s; }
    std::vector<std::string>& getSpectrum() { return spectrum; }
    const std::vector<std::string>& getSpectrum() const { return spectrum; }
    int findStartVertexIndex(const DNAInstance& instance);


    
};

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
    
    // Add missing method declarations
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
 * Klasa do zapisu i odczytu instancji z/do pliku.
 */
class InstanceIO
{
public:
    static bool saveInstance(const DNAInstance &instance, const std::string &filename);
    static bool loadInstance(const std::string &filename, DNAInstance &instance);
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
 * Budowniczy tworzący obiekt DNAInstance:
 *  - ustala parametry (n, k, deltaK, lNeg, lPoz, repAllowed, probablePositive)
 *  - generuje DNA i spektrum
 *  - nakłada błędy (jeśli lNeg > 0 lub lPoz > 0)
 */
class DNAInstanceBuilder {
private:
    DNAInstance instance;    // budowany obiekt
    ErrorContext errorCtx;   // kontekst do nakładania strategii błędów

public:
    // Settery łańcuchowe parametrów
    DNAInstanceBuilder& setN(int n);
    DNAInstanceBuilder& setK(int k);
    DNAInstanceBuilder& setDeltaK(int dk);
    DNAInstanceBuilder& setLNeg(int ln);
    DNAInstanceBuilder& setLPoz(int lp);
    DNAInstanceBuilder& setRepAllowed(bool rep);
    DNAInstanceBuilder& setProbablePositive(int val);

    // Budowanie DNA i spektrum
    DNAInstanceBuilder& buildDNA();
    DNAInstanceBuilder& buildSpectrum();

    // Nakładanie błędów
    DNAInstanceBuilder& applyError(IErrorIntroductionStrategy* strategy);

    // Zwrócenie gotowej instancji
    DNAInstance getInstance();
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

#endif //DNA_GENERATOR_H
