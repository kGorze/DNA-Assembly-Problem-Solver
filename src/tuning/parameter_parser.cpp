//
// Created by konrad_guest on 23/01/2025.
//
#include "../include/tuning/parameters_parser.h"
#include <algorithm>
#include <random>

// Implementacja parseConfigFile
ParameterSet ParameterParser::parseConfigFile(const std::string &filename) {
    ParameterSet ps;
    std::ifstream fin(filename);
    if (!fin.is_open()) {
        std::cerr << "Could not open param config file: " << filename << "\n";
        return ps;
    }
    std::string line;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        auto pos = line.find(',');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string val = line.substr(pos + 1);
            ps.params[key] = val;
        }
    }
    return ps;
}

// Implementacja generateRandomCandidates
std::vector<ParameterSet> ParameterParser::generateRandomCandidates(int numCandidates) {
    std::vector<ParameterSet> result;
    result.reserve(numCandidates);

    // Te same zakresy co wyżej (lub węższe):
    std::vector<int> popSizes = {20,50,80,100,130,150,200};
    std::vector<double> replacementRatios = {0.1,0.2,0.3,0.4,0.5,0.6};
    std::vector<int> tournamentSizes = {3,5,10};
    std::vector<double> mutRates = {0.05, 0.1, 0.2, 0.3,0.4};

    std::vector<double> inertias = {0.5, 0.6, 0.7, 0.8};
    std::vector<int> intervals = {5, 10, 20};
    std::vector<int> minTrials = {3, 5, 10};
    std::vector<double> minProbs = {0.05, 0.1, 0.2};

    // Silnik losowy
    std::mt19937 rng(std::random_device{}());

    // Rozkłady
    std::uniform_int_distribution<int> distPop(0, (int)popSizes.size()-1);
    std::uniform_int_distribution<int> distRR(0, (int)replacementRatios.size()-1);
    std::uniform_int_distribution<int> distTS(0, (int)tournamentSizes.size()-1);
    std::uniform_int_distribution<int> distMut(0, (int)mutRates.size()-1);
    std::uniform_int_distribution<int> distInertia(0, (int)inertias.size()-1);
    std::uniform_int_distribution<int> distInterval(0, (int)intervals.size()-1);
    std::uniform_int_distribution<int> distMinTrials(0, (int)minTrials.size()-1);
    std::uniform_int_distribution<int> distMinProbs(0, (int)minProbs.size()-1);

    for(int i = 0; i < numCandidates; i++) {
        ParameterSet pset;
        pset.params["populationSize"] = std::to_string(popSizes[distPop(rng)]);
        pset.params["mutationRate"]   = std::to_string(mutRates[distMut(rng)]);
        pset.params["replacementRatio"] = std::to_string(replacementRatios[distRR(rng)]);
        pset.params["tournamentSize"]  = std::to_string(tournamentSizes[distTS(rng)]);

        pset.params["crossoverType"] = "adaptive";
        pset.params["adaptive.inertia"] = std::to_string(inertias[distInertia(rng)]);
        pset.params["adaptive.adaptationInterval"] = std::to_string(intervals[distInterval(rng)]);
        pset.params["adaptive.minTrials"] = std::to_string(minTrials[distMinTrials(rng)]);
        pset.params["adaptive.minProb"]   = std::to_string(minProbs[distMinProbs(rng)]);

        pset.params["selectionMethod"] = "tournament";

        result.push_back(pset);
    }

    return result;
}

// Implementacja generateGridOfCandidates
std::vector<ParameterSet> ParameterParser::generateGridOfCandidatesWithout() {
    std::vector<ParameterSet> candidates;
    
    // Define parameter ranges
    std::vector<int> populationSizes = {50, 100, 200};
    std::vector<int> maxGenerations = {100, 500, 1000};
    std::vector<double> mutationRates = {0.01, 0.05, 0.1, 0.2};
    std::vector<double> crossoverRates = {0.6, 0.7, 0.8, 0.9};
    std::vector<int> tournamentSizes = {2, 3, 4, 5};
    
    // Generate all combinations
    for (int popSize : populationSizes) {
        for (int maxGen : maxGenerations) {
            for (double mutRate : mutationRates) {
                for (double crossRate : crossoverRates) {
                    for (int tournSize : tournamentSizes) {
                        ParameterSet params;
                        params.populationSize = popSize;
                        params.maxGenerations = maxGen;
                        params.mutationRate = mutRate;
                        params.crossoverRate = crossRate;
                        params.tournamentSize = tournSize;
                        candidates.push_back(params);
                    }
                }
            }
        }
    }
    
    // Shuffle candidates to avoid bias
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(candidates.begin(), candidates.end(), gen);
    
    return candidates;
}

// Implementacja generateGridOfCandidatesRS
std::vector<ParameterSet> ParameterParser::generateGridOfCandidatesRS() {
    int numCandidates = 1000000; // Ustal liczbę kandydatów
    return generateRandomCandidates(numCandidates);
}