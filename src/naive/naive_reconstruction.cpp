//
// Created by konrad_guest on 29/12/2024.
//

#include "naive/naive_reconstruction.h"
#include <algorithm>  // std::sort, std::shuffle
#include <random>
#include <sstream>
#include <chrono>  

std::string NaiveReconstructor::reconstructDNA(const DNAInstance &instance,
                                               NaiveReconstructionMethod method)
{
    switch(method) {
        case NaiveReconstructionMethod::METHOD_A:
            return reconstructDNA_A(instance);
        case NaiveReconstructionMethod::METHOD_B:
            return reconstructDNA_B(instance);
        case NaiveReconstructionMethod::METHOD_C:
            return reconstructDNA_C(instance);
        default:
            // Domyślnie (lub dla błędnej wartości) - np. metoda A
            return reconstructDNA_A(instance);
    }
}

// ================== Metoda A ==================
std::string NaiveReconstructor::reconstructDNA_A(const DNAInstance &instance)
{
    // 1. Kopiujemy spektrum
    std::vector<std::string> spectrum = instance.getSpectrum();
    int k = instance.getK();
    if(spectrum.empty()) {
        return "";
    }

    // 2. Sortowanie alfabetyczne
    std::sort(spectrum.begin(), spectrum.end());

    // 3. Naiwne sklejanie (sort + doklejanie ostatniego znaku)
    std::ostringstream oss;
    oss << spectrum[0];  // pierwszy k-mer w całości
    for (size_t i = 1; i < spectrum.size(); ++i) {
        // Doklejamy ostatni znak (lub kilka znaków)
        oss << spectrum[i].substr(k-1);
    }

    return oss.str();
}

// ================== Metoda B ==================
std::string NaiveReconstructor::reconstructDNA_B(const DNAInstance &instance)
{
    std::vector<std::string> spectrum = instance.getSpectrum();
    int k = instance.getK();
    if(spectrum.empty()) {
        return "";
    }

    // Weźmy pierwszy k-mer
    std::string result = spectrum[0];
    
    // Iterujemy przez resztę k-merów w ich naturalnej kolejności
    for (size_t i = 1; i < spectrum.size(); ++i) {
        const std::string &current = spectrum[i];
        // Znajdź największe dopasowanie sufiks->prefiks
        int max_overlap = 0;
        for(int overlap = k-1; overlap > 0; --overlap) {
            if(result.compare(result.size() - overlap, overlap, current, 0, overlap) == 0) {
                max_overlap = overlap;
                break;
            }
        }
        // Dopisz resztę `current` do `result`
        result += current.substr(max_overlap);
    }

    return result;
}

// ================== Metoda C ==================
// Dla przykładu: losujemy kolejność k-merów i sklejamy
// w podobny sposób jak w A (doklejanie kilku znaków).
std::string NaiveReconstructor::reconstructDNA_C(const DNAInstance &instance)
{
    std::vector<std::string> spectrum = instance.getSpectrum();
    int k = instance.getK();
    if(spectrum.empty()) {
        return "";
    }

    // Wylosuj kolejność k-merów
    auto &rng = RandomGenerator::getInstance().get();
    std::shuffle(spectrum.begin(), spectrum.end(), rng);

    // Sklej tak jak w metodzie A, tylko w losowej kolejności
    std::ostringstream oss;
    oss << spectrum[0];
    for (size_t i = 1; i < spectrum.size(); ++i) {
        oss << spectrum[i].substr(k-1);
    }
    return oss.str();
}

int levenshteinDistance(const std::string &s1, const std::string &s2) {
    int len1 = (int)s1.size();
    int len2 = (int)s2.size();
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1, 0));

    for(int i = 0; i <= len1; ++i) dp[i][0] = i;
    for(int j = 0; j <= len2; ++j) dp[0][j] = j;

    for(int i = 1; i <= len1; ++i) {
        for(int j = 1; j <= len2; ++j) {
            int cost = (s1[i-1] == s2[j-1]) ? 0 : 1;
            dp[i][j] = std::min({ 
                dp[i-1][j] + 1, 
                dp[i][j-1] + 1, 
                dp[i-1][j-1] + cost
            });
        }
    }
    return dp[len1][len2];
}