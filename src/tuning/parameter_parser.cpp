//
// Created by konrad_guest on 23/01/2025.
//
#include "../include/tuning/parameters_parser.h"

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
    std::vector<ParameterSet> result;
    
    // Reduce parameter ranges to prevent memory exhaustion
    std::vector<int> popSizes = {50, 100, 200, 300};
    std::vector<double> replacementRatios = {0.2, 0.4, 0.6};
    std::vector<int> tournamentSizes = {3, 5, 10};
    std::vector<double> mutRates = {0.1, 0.3, 0.5};
    std::vector<double> inertias = {0.6, 0.7};
    std::vector<int> intervals = {5, 10};
    std::vector<int> minTrials = {3, 5};
    std::vector<double> minProbs = {0.1, 0.2};

    // Pre-calculate total size and reserve memory
    size_t totalCombinations = popSizes.size() * replacementRatios.size() * 
                              tournamentSizes.size() * mutRates.size() * 
                              inertias.size() * intervals.size() * 
                              minTrials.size() * minProbs.size();
    
    // Reserve space to prevent reallocation
    result.reserve(totalCombinations);

    // Generate combinations with proper initialization
    for (int ps : popSizes) {
        for (double mr : mutRates) {
            for (double rr : replacementRatios) {
                for (int ts : tournamentSizes) {
                    for (double in : inertias) {
                        for (int ai : intervals) {
                            for (int mt : minTrials) {
                                for (double mp : minProbs) {
                                    ParameterSet pset;
                                    
                                    // Initialize map with exact size needed
                                    pset.params.reserve(10);
                                    
                                    // Use emplace to avoid extra copies
                                    pset.params.emplace("populationSize", std::to_string(ps));
                                    pset.params.emplace("mutationRate", std::to_string(mr));
                                    pset.params.emplace("replacementRatio", std::to_string(rr));
                                    pset.params.emplace("tournamentSize", std::to_string(ts));
                                    pset.params.emplace("crossoverType", "adaptive");
                                    pset.params.emplace("adaptive.inertia", std::to_string(in));
                                    pset.params.emplace("adaptive.adaptationInterval", std::to_string(ai));
                                    pset.params.emplace("adaptive.minTrials", std::to_string(mt));
                                    pset.params.emplace("adaptive.minProb", std::to_string(mp));
                                    pset.params.emplace("selectionMethod", "tournament");
                                    
                                    result.push_back(std::move(pset));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    return result;
}

// Implementacja generateGridOfCandidatesRS
std::vector<ParameterSet> ParameterParser::generateGridOfCandidatesRS() {
    int numCandidates = 1000000; // Ustal liczbę kandydatów
    return generateRandomCandidates(numCandidates);
}