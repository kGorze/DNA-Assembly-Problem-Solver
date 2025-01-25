//
// Created by konrad_guest on 23/01/2025.
//

#ifndef TUNING_STRUCTURES_H
#define TUNING_STRUCTURES_H

#include <string>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <stdexcept>

/**
 * Reprezentacja zestawu parametrów do tuningu.
 * Może zawierać np. parametry liczbowe, ciągłe, symboliczne.
 * Przykładowo: "populationSize=100", "mutationRate=0.1", "selectionMethod=tournament" itd.
 */
struct ParameterSet {
    // Aby ułatwić, możemy przechowywać w mapie string->string,
    // a interpretacją typów (int, double, itp.) zająć się osobno.
    std::unordered_map<std::string, std::string> params;
    
    bool contains(const std::string& key) const {
        return params.find(key) != params.end();
    }
    
    int getInt(const std::string& key) const {
        auto it = params.find(key);
        if (it != params.end()) {
            return std::stoi(it->second);
        }
        throw std::runtime_error("Parameter not found: " + key);
    }
    
    double getDouble(const std::string& key) const {
        auto it = params.find(key);
        if (it != params.end()) {
            return std::stod(it->second);
        }
        throw std::runtime_error("Parameter not found: " + key);
    }
    
    std::string getString(const std::string& key) const {
        auto it = params.find(key);
        if (it != params.end()) {
            return it->second;
        }
        throw std::runtime_error("Parameter not found: " + key);
    }
    
    // Funkcja pomocnicza do wypisywania w CSV, debug, itp.
    std::string toString() const {
        std::string result;
        for (auto &p : params) {
            result += p.first + "=" + p.second + ";";
        }
        return result;
    }
};

/**
 * Struktura wyników pojedynczego uruchomienia AE z zadanym zestawem parametrów.
 */
struct TuningResult {

    
    ParameterSet parameterSet;
    double fitness = 0.0;       // Ostateczny fitness
    double executionTime = 0.0; // Czas wykonania (sekundy albo ms)
    
    // Można dodać dowolne metryki – np. coverage, edgeScore, itp.
    // "Dodatkowe metryki wydajności":
    std::unordered_map<std::string, double> extraMetrics;
    
    std::string toCSV() const {
        // Przykład – finalny format można zmodyfikować
        // Kolumny: "Testowane parametry; Uzyskany fitness; Czas wykonania; [Extra1]; [Extra2]..."
        // Dla łatwego wczytania w Excelu możemy użyć np. ';' jako separatora
        std::string csvLine;
        csvLine += "\"" + parameterSet.toString() + "\";";
        csvLine += std::to_string(fitness) + ";";
        csvLine += std::to_string(executionTime);
        // Dopisanie extra metryk
        for (auto &m : extraMetrics) {
            csvLine += ";" + m.first + "=" + std::to_string(m.second);
        }
        return csvLine;
    }
};



#endif //TUNING_STRUCTURES_H
