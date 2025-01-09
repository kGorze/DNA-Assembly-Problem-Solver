//
// Created by konrad_guest on 07/01/2025.
//
#include <random>
#include <algorithm>
#include <iostream>
#include "metaheuristics/representation.h"

//TODO
// - W pełnym kodzie należałoby jeszcze zadbać o zwalnianie pamięci (np. w destruktorze algorytmu) oraz bardziej wyrafinowane inicjalizacje i dekodowania. Powyższe pokazuje ideę.

/************************************************
 *  DirectDNARepresentation
 ************************************************/
std::vector<void*> DirectDNARepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<void*> population;
    population.reserve(popSize);

    // Uzyskaj rzeczywistą długość n z instancji
    int n = instance.getN();
    if(n <= 0) n = fallbackN;  // fallback w przypadku problemów

    int k = instance.getK();
    const auto &spectrum = instance.getSpectrum();
    if(spectrum.empty() || k <= 0) {
        // Fallback: gdy spektrum jest puste, generuj losowe DNA
        // (możesz tutaj użyć wcześniejszego kodu)
        return population;
    }

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<size_t> spectrumDist(0, spectrum.size() - 1);
    static const char* letters = "ACGT";

    for(int i = 0; i < popSize; i++) {
        auto *dnaStr = new std::string;
        dnaStr->reserve(n);

        // 2.1. Losowo wybierz pierwszy k-mer ze spektrum
        std::string current = spectrum[spectrumDist(rng)];
        *dnaStr = current;

        // 2.2. Łączenie k-merów aż do osiągnięcia długości n
        while((int)dnaStr->size() < n) {
            int maxOverlap = 0;
            std::string bestNext;
            // Szukaj k-mera z największym dopasowaniem sufiks->prefiks
            for(const auto &km : spectrum) {
                // Ograniczenie długości nakładki do min(k-1, aktualna długość dnaStr)
                int maxPossibleOverlap = std::min(k - 1, (int)dnaStr->size());
                for(int overlap = maxPossibleOverlap; overlap > 0; overlap--) {
                    // Porównanie sufiksu current dnaStr z prefiksem km
                    if(dnaStr->compare(dnaStr->size() - overlap, overlap, km, 0, overlap) == 0) {
                        if(overlap > maxOverlap) {
                            maxOverlap = overlap;
                            bestNext = km;
                        }
                        break; // znaleziono dopasowanie dla tego km, przejdź do następnego
                    }
                }
            }
            if(maxOverlap == 0) {
                // Jeżeli nie znaleziono żadnego k-mera pasującego do sufiksu,
                // dodaj losowy znak (alternatywnie: wstrzymaj lub zakończ)
                dnaStr->push_back(letters[rng() % 4]);
            } else {
                // Dopisz resztę k-mera, która się nie pokrywa
                dnaStr->append(bestNext.substr(maxOverlap));
            }
        }

        population.push_back(dnaStr);
    }
    return population;
}

std::string DirectDNARepresentation::decodeToDNA(void* individual, const DNAInstance &instance)
{
    // Tutaj "individual" to w rzeczywistości (std::string*)
    auto dnaPtr = static_cast<std::string*>(individual);
    return *dnaPtr; 
}

/************************************************
 *  PermutationRepresentation
 ************************************************/
std::vector<void*> PermutationRepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<void*> population;
    population.reserve(popSize);

    // W tej reprezentacji osobnik to permutacja indeksów k-merów, np.
    // (0,1,2,3,4,...) w pewnej kolejności, bądź same k-mery.
    // Na potrzeby przykładu zrobimy permutację [0..99].
    // (Realnie – w SBH: permutację k-merów z instance.getSpectrum().)
    int length = instance.getSpectrum().size();  
    for(int i = 0; i < popSize; i++){
        auto *perm = new std::vector<int>;
        perm->reserve(length);
        for(int idx = 0; idx < length; idx++){
            perm->push_back(idx);
        }
        // Wylosuj kolejność
        std::shuffle(perm->begin(), perm->end(), std::mt19937(std::random_device{}()));
        population.push_back(perm);
    }
    return population;
}

std::string PermutationRepresentation::decodeToDNA(void* individual, const DNAInstance &instance)
{
    if (individual == nullptr) {
        return "";
    }

    auto permPtr = static_cast<std::vector<int>*>(individual);
    if (!permPtr || permPtr->empty()) {
        return "";
    }

    const auto &spec = instance.getSpectrum();
    if (spec.empty()) {
        return "";
    }

    int k = instance.getK();
    if (k <= 0) {
        return "";
    }

    int firstIndex = (*permPtr)[0];
    if (firstIndex < 0 || firstIndex >= static_cast<int>(spec.size())) {
        return "";
    }

    std::string result = spec[firstIndex];

    for (size_t i = 1; i < permPtr->size(); ++i) {
        int idx = (*permPtr)[i];
        if (idx < 0 || idx >= static_cast<int>(spec.size())) {
            continue;
        }

        const std::string &km = spec[idx];
        if (km.size() < static_cast<size_t>(k-1)) {
            continue;
        }

        result += km.substr(k-1);
    }

    return result;
}

/************************************************
 *  GraphPathRepresentation
 ************************************************/
std::vector<void*> GraphPathRepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<void*> population;
    population.reserve(popSize);
    // Uproszczone: "ścieżka w grafie" będzie wektorem int
    // (podobnie do PermutationRepresentation, ale z inną interpretacją).
    int length = 50;  
    for(int i = 0; i < popSize; i++){
        auto *path = new std::vector<int>;
        path->reserve(length);
        for(int idx = 0; idx < length; idx++){
            path->push_back(idx);
        }
        // Możemy wylosować kolejność lub wylosować dokąd prowadzą krawędzie.
        std::shuffle(path->begin(), path->end(), std::mt19937(std::random_device{}()));
        population.push_back(path);
    }
    return population;
}

std::string GraphPathRepresentation::decodeToDNA(void* individual, const DNAInstance &instance)
{
    // W praktyce używamy "overlap graph" lub "de Bruijn graph".
    // Uproszczone – podobnie do permutacji, tylko interpretacja inna
    // (zakładamy, że "vector<int>" to ścieżka odwiedzenia węzłów w grafie).
    auto pathPtr = static_cast<std::vector<int>*>(individual);
    const auto &spec = instance.getSpectrum();
    if(spec.empty() || pathPtr->empty()) {
        return "";
    }
    // Bierzemy pierwszy "węzeł" i doklejamy dalej
    // Kod identyczny do PermutationRepresentation, ale logicznie to "ścieżka w grafie".
    std::string result = spec[(*pathPtr)[0]];
    int k = instance.getK();
    for(size_t i = 1; i < pathPtr->size(); i++){
        int idx = (*pathPtr)[i];
        if(idx < 0 || idx >= (int)spec.size()) {
            continue;
        }
        result += spec[idx].substr(k-1);
    }
    return result;
}