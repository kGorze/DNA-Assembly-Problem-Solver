//
// Created by konrad_guest on 07/01/2025.
// SMART

#include "metaheuristics/representation.h"
#include <random>
#include <algorithm>
#include <iostream>

/************************************************
 *  DirectDNARepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
DirectDNARepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    int n = instance.getN();
    if (n <= 0) n = fallbackN; 

    std::mt19937 rng(std::random_device{}());
    static const char letters[] = "ACGT";

    for (int i = 0; i < popSize; i++) {
        auto individual = std::make_shared<std::vector<int>>();
        individual->reserve(n);

        for (int j = 0; j < n; j++) {
            int base = rng() % 4; // 0..3
            individual->push_back(base);
        }
        population.push_back(individual);
    }

    return population;
}

std::string 
DirectDNARepresentation::decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                                     const DNAInstance &instance)
{
    static const char letters[] = {'A','C','G','T'};
    if (!individual) return "";

    std::string dna;
    dna.reserve(individual->size());
    for (auto val : *individual) {
        int idx = std::max(0, std::min(val, 3));
        dna.push_back(letters[idx]);
    }

    return dna;
}

/************************************************
 *  PermutationRepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
PermutationRepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    const auto &spec = instance.getSpectrum();
    int length = (int)spec.size();
    if (length <= 0) {
        return population;
    }

    for(int i = 0; i < popSize; i++){
        auto perm = std::make_shared<std::vector<int>>();
        perm->reserve(length);
        for(int idx = 0; idx < length; idx++){
            perm->push_back(idx);
        }
        std::shuffle(perm->begin(), perm->end(), std::mt19937(std::random_device{}()));
        population.push_back(perm);
    }
    return population;
}

std::string 
PermutationRepresentation::decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                                       const DNAInstance &instance)
{
    if (!individual || individual->empty()) {
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

    // Sprawdzamy, czy rozmiar permutacji zgadza się z rozmiarem spektrum
    if ((int)individual->size() != (int)spec.size()) {
        // np. ostrzeżenie
        //std::cerr << "[PermutationRepresentation] Warning: mismatch sizes.\n";
        // i wczytujemy minimalnie
    }

    int firstIndex = (*individual)[0];
    if (firstIndex < 0 || firstIndex >= (int)spec.size()) {
        return ""; // błąd
    }

    std::string result = spec[firstIndex];

    for (size_t i = 1; i < individual->size(); i++) {
        int idx = (*individual)[i];
        if (idx < 0 || idx >= (int)spec.size()) {
            // pomijamy
            continue;
        }
        const std::string &km = spec[idx];
        if ((int)km.size() < k - 1) {
            continue;
        }
        // doklejamy last (k-1) znaków
        result += km.substr(k - 1);
    }

    return result;
}

/************************************************
 *  GraphPathRepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
GraphPathRepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    int length = 50;  // w oryg. stała – do poprawy dla lepszej jakości
    for(int i = 0; i < popSize; i++) {
        auto path = std::make_shared<std::vector<int>>();
        path->reserve(length);
        for(int idx = 0; idx < length; idx++) {
            path->push_back(idx);
        }
        std::shuffle(path->begin(), path->end(), std::mt19937(std::random_device{}()));
        population.push_back(path);
    }

    return population;
}

std::string 
GraphPathRepresentation::decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                                     const DNAInstance &instance)
{
    if (!individual || individual->empty()) {
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

    // Podobnie jak w PermutationRepresentation
    int firstIndex = (*individual)[0];
    if (firstIndex < 0 || firstIndex >= (int)spec.size()) {
        return "";
    }

    std::string result = spec[firstIndex];

    for (size_t i = 1; i < individual->size(); i++) {
        int idx = (*individual)[i];
        if (idx < 0 || idx >= (int)spec.size()) {
            continue;
        }
        result += spec[idx].substr(k - 1);
    }

    return result;
}
