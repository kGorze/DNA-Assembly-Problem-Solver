//
// Created by konrad_guest on 28/12/2024.
//

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


class DNAInstance {
private:
    std::string dna;
    int n;
    int k;
    std::vector<std::string> spectrum;
public:
    // get/set
    void setDNA(const std::string &d) { dna = d; }
    std::string getDNA() const { return dna; }

    void setN(int val) { n = val; }
    int getN() const { return n; }

    void setK(int val) { k = val; }
    int getK() const { return k; }

    void setSpectrum(const std::vector<std::string> &s) { spectrum = s; }
    std::vector<std::string>& getSpectrum() { return spectrum; }
    const std::vector<std::string>& getSpectrum() const { return spectrum; }
};

class DNAGenerator {
public:
    /**
     * Generuje łańcuch DNA o długości n, zwracając go w postaci std::string.
     */
    std::string generateDNA(int n);
};

/**
 * Klasa odpowiedzialna za generowanie spektrum (k-merów)
 * na podstawie przekazanego łańcucha DNA i parametru k.
 */
class SpectrumGenerator {
public:
    /**
     * Zwraca wektor wszystkich k-merów (oligonukleotydów) o długości k
     * dla podanego łańcucha DNA.
     */
    std::vector<std::string> generateSpectrum(const std::string &dna, int k);
};


class InstanceIO
{
public:
    // Statyczna metoda zapisu, zwraca true w razie sukcesu, false w przypadku błędu
    static bool saveInstance(const DNAInstance &instance, const std::string &filename);

    // Statyczna metoda odczytu, zwraca true w razie sukcesu, false w przypadku błędu
    // Odczytane dane zapisuje w 'instance'
    static bool loadInstance(const std::string &filename, DNAInstance &instance);
};

class RandomGenerator {
private:
    std::mt19937 gen;
    static RandomGenerator* instance;

    RandomGenerator() {
        // np. seed oparty o std::random_device + high_resolution_clock
        std::random_device rd;
        gen.seed(rd() ^ (unsigned long long)std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }

public:
    static RandomGenerator& getInstance() {
        static RandomGenerator instance;
        return instance;
    }

    std::mt19937& get() {
        return gen;
    }

    // usuń kopiowanie, aby klasa faktycznie była Singletonem
    RandomGenerator(const RandomGenerator&) = delete;
    RandomGenerator& operator=(const RandomGenerator&) = delete;
};

class IErrorIntroductionStrategy {
public:
    virtual ~IErrorIntroductionStrategy() = default;

    // Funkcja wprowadzająca błędy do podanego obiektu DNAInstance.
    // Konkretne implementacje (np. NegativeErrorIntroducer, PositiveErrorIntroducer)
    // będą nadpisywać tę metodę i modyfikować spektrum według własnej logiki.
    virtual void introduceErrors(DNAInstance &instance) = 0;
};

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


class DNAInstanceBuilder {
private:
    DNAInstance instance;    // budowany obiekt
    ErrorContext errorCtx;   // do nakładania strategii błędów

public:
    DNAInstanceBuilder& setN(int n) {
        instance.setN(n);
        return *this;
    }

    DNAInstanceBuilder& setK(int k) {
        instance.setK(k);
        return *this;
    }

    DNAInstanceBuilder& buildDNA() {
        // Wygeneruj DNA i zapisz do instance
        DNAGenerator dnaGen;
        std::string dna = dnaGen.generateDNA(instance.getN());
        instance.setDNA(dna);
        return *this;
    }

    DNAInstanceBuilder& buildSpectrum() {
        SpectrumGenerator specGen;
        auto spec = specGen.generateSpectrum(instance.getDNA(), instance.getK());
        instance.setSpectrum(spec);
        return *this;
    }

    DNAInstanceBuilder& applyError(IErrorIntroductionStrategy* strategy) {
        // Deleguje wprowadzanie błędów do obiektu ErrorContext
        errorCtx.setStrategy(strategy);
        errorCtx.execute(instance);
        return *this;
    }

    DNAInstance getInstance() {
        return instance;
    }
};


class NegativeErrorIntroducer : public IErrorIntroductionStrategy {
private:
    double percentNeg;
public:
    NegativeErrorIntroducer(double p): percentNeg(p) {}

    virtual void introduceErrors(DNAInstance &instance) override {
        // Usuwamy z wektora 'spectrum' pewien % elementów
        auto &spectrum = instance.getSpectrum(); // referencja
        // Ile usunąć
        int toRemove = (int)(spectrum.size() * (percentNeg / 100.0) + 0.5);

        // Wylosuj kolejność
        auto &rng = RandomGenerator::getInstance().get();
        std::shuffle(spectrum.begin(), spectrum.end(), rng);

        if(toRemove > (int)spectrum.size()) {
            toRemove = (int)spectrum.size();
        }
        spectrum.erase(spectrum.begin(), spectrum.begin() + toRemove);
    }
};

class PositiveErrorIntroducer : public IErrorIntroductionStrategy {
private:
    double percentPos;
    double newKmerProbability; // Prawdopodobieństwo dodania nowego k-meru

public:
    // Konstruktor z dodatkowym parametrem określającym prawdopodobieństwo dodania nowego k-meru
    PositiveErrorIntroducer(double p, double newProb = 0.7) 
        : percentPos(p), newKmerProbability(newProb) {}

    virtual void introduceErrors(DNAInstance &instance) override {
        auto &spectrum = instance.getSpectrum();
        int total = static_cast<int>(spectrum.size());
        int toAdd = static_cast<int>(total * (percentPos / 100.0) + 0.5);

        static const std::string nucleotides = "ACGT";
        auto &rng = RandomGenerator::getInstance().get();
        std::uniform_real_distribution<double> probDist(0.0, 1.0);
        std::uniform_int_distribution<int> nucDist(0, 3);

        // Zbiór istniejących fragmentów
        std::unordered_set<std::string> existing(spectrum.begin(), spectrum.end());

        for(int i = 0; i < toAdd; ++i) {
            double prob = probDist(rng);
            std::string kmer;

            if(prob < newKmerProbability && static_cast<int>(existing.size()) < pow(4, instance.getK())) {
                // Dodajemy nowy k-mer
                while(true) {
                    kmer.clear();
                    for(int j = 0; j < instance.getK(); ++j) {
                        kmer += nucleotides[nucDist(rng)];
                    }
                    if(existing.find(kmer) == existing.end()) {
                        spectrum.push_back(kmer);
                        existing.insert(kmer);
                        break;
                    }
                    // Jeśli wygenerowany k-mer już istnieje, próbujemy ponownie
                }
            } else {
                // Dodajemy duplikat istniejącego k-meru
                if(!spectrum.empty()) {
                    std::uniform_int_distribution<int> indexDist(0, static_cast<int>(spectrum.size()) - 1);
                    int randomIndex = indexDist(rng);
                    kmer = spectrum[randomIndex];
                    spectrum.push_back(kmer);
                }
            }
        }
    }
};



#endif //DNA_GENERATOR_H
