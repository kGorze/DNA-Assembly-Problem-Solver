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

    // PRZYKŁADOWA inicjalizacja: generujemy losowy łańcuch DNA stałej długości,
    // np. n=500. W practice – powinniśmy znać n skądś, np. z parametru lub z DNAInstance.
    // Dla uproszczenia tu napisane "500".
    int n = 500; 
    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<int> dist(0, 3);
    static const char* letters = "ACGT";

    for(int i = 0; i < popSize; i++) {
        // Tworzymy nowy std::string
        auto *dnaStr = new std::string;
        dnaStr->reserve(n);
        for(int j = 0; j < n; j++){
            (*dnaStr) += letters[dist(rng)];
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
    int length = 100;  
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
    // Tutaj "individual" to (std::vector<int>*). 
    // Dekodowanie: 
    //  1. Weź spektrum k-merów z instance.getSpectrum().
    //  2. Iteruj po kolei po indeksach w permutacji.
    //  3. "Sklejaj" k-mery z maksymalnym overlapem (tu w dużym uproszczeniu).
    auto permPtr = static_cast<std::vector<int>*>(individual);
    const auto &spec = instance.getSpectrum();
    if(spec.empty() || permPtr->empty()) {
        return "";
    }
    // Bierzemy pierwszy k-mer
    std::string result = spec[ (*permPtr)[0] ];
    int k = instance.getK();
    for(size_t i = 1; i < permPtr->size(); i++){
        int idx = (*permPtr)[i];
        if(idx < 0 || idx >= (int)spec.size()) {
            // W praktyce: błąd w rejonach permutacji
            continue;
        }
        // Doklej overlap
        // (naiwne: doklejamy ostatnie k-1 znaków)
        result += spec[idx].substr(k-1);
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