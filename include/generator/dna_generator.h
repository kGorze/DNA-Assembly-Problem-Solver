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

    std::vector<std::string> spectrum;

public:
    std::vector<double> genes;  // Keep the existing implementation

    
    // get/set
    void setDNA(const std::string &d) { dna = d; }
    std::string getDNA() const { return dna; }

    void setN(int val) { n = val; }
    int getN() const { return n; }

    void setK(int val) { k = val; }
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
};

/**
 * Klasa generująca DNA (ciąg znaków 'ACGT') o zadanej długości n.
 * Jeżeli repAllowed = false, można tu rozszerzyć logikę ograniczającą
 * częste powtórzenia (w przykładzie pomijamy, bo zależy mocno od definicji
 * "powtórzeń").
 */
class DNAGenerator {
public:
    std::string generateDNA(int n, bool repAllowed);
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
