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
#include <initializer_list>

/**
 * Reprezentacja zestawu parametrów do tuningu.
 * Może zawierać np. parametry liczbowe, ciągłe, symboliczne.
 * Przykładowo: "populationSize=100", "mutationRate=0.1", "selectionMethod=tournament" itd.
 */
struct ParameterSet {
    std::unordered_map<std::string, std::string> params;
    
    // Add default constructor
    ParameterSet() {
        // Reserve space for typical number of parameters
        params.reserve(10);
    }
    
    // Add copy constructor
    ParameterSet(const ParameterSet& other) {
        params.reserve(other.params.size());
        params = other.params;
    }
    
    // Add move constructor
    ParameterSet(ParameterSet&& other) noexcept 
        : params(std::move(other.params)) {}
    
    // Add assignment operator
    ParameterSet& operator=(const ParameterSet& other) {
        if (this != &other) {
            params.clear();
            params.reserve(other.params.size());
            params = other.params;
        }
        return *this;
    }
    
    // Add move assignment operator
    ParameterSet& operator=(ParameterSet&& other) noexcept {
        if (this != &other) {
            params = std::move(other.params);
        }
        return *this;
    }
    
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

    // Getters for common parameters with default values
    int getPopulationSize() const { return contains("populationSize") ? getInt("populationSize") : 100; }
    int getMaxGenerations() const { return contains("maxGenerations") ? getInt("maxGenerations") : 100; }
    double getMutationRate() const { return contains("mutationRate") ? getDouble("mutationRate") : 0.1; }
    double getCrossoverRate() const { return contains("crossoverRate") ? getDouble("crossoverRate") : 0.8; }
    int getTournamentSize() const { return contains("tournamentSize") ? getInt("tournamentSize") : 5; }

    // Setters for common parameters
    void setPopulationSize(int value) { params["populationSize"] = std::to_string(value); }
    void setMaxGenerations(int value) { params["maxGenerations"] = std::to_string(value); }
    void setMutationRate(double value) { params["mutationRate"] = std::to_string(value); }
    void setCrossoverRate(double value) { params["crossoverRate"] = std::to_string(value); }
    void setTournamentSize(int value) { params["tournamentSize"] = std::to_string(value); }
    
    // Funkcja pomocnicza do wypisywania w CSV, debug, itp.
    std::string toString() const {
        std::string result;
        for (auto &p : params) {
            result += p.first + "=" + p.second + ";";
        }
        return result;
    }
    
    // Add destructor
    ~ParameterSet() = default;
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

    // Constructor for easy initialization
    TuningResult() = default;
    
    TuningResult(const ParameterSet& ps, double fit, double time)
        : parameterSet(ps)
        , fitness(fit)
        , executionTime(time)
        , extraMetrics() {}
    
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
