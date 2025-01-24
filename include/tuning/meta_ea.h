//
// Created by konrad_guest on 23/01/2025.
//

#ifndef META_EA_H
#define META_EA_H

#include "racing.h"  // Dodaj ten nagłówek
#include "../generator/dna_generator.h"
#include "configuration/genetic_algorithm_configuration.h" // Dołączenie GAConfig
#include "metaheuristics/genetic_algorithm.h"  

#include "tuning_structures.h"
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>

struct MetaEAConfig {
    int populationSize = 20;
    int maxGenerations = 10;
    void validate() const {
        if (populationSize < 2) 
            throw std::invalid_argument("Population size must be at least 2");
        if (maxGenerations < 1)
            throw std::invalid_argument("Generations must be at least 1");
    }
};

/**
 * Przykładowa klasa MetaEA, która optymalizuje parametry AE.
 */
class MetaEA {
public:
    MetaEA(const MetaEAConfig &cfg, const Racing::Configuration &rcfg) // Korzystanie z namespace
       : config(cfg), racingCfg(rcfg) {
        config.validate(); // Validate configuration
    }

    /**
     * Funkcja uruchamia Meta-EA, aby znaleźć najlepsze parametry AE.
     * Wewnątrz, w każdej generacji, możemy wykorzystać Racing
     * lub normalną ocenę (wielokrotne uruchomienia) – w zależności od koncepcji.
     */
    ParameterSet runMetaEA(std::function<TuningResult(const ParameterSet&)> evaluator)
    {
        // Get singleton instance
        auto& config = GAConfig::getInstance();
        
        try {
            config.loadFromFile("config.cfg");
        } catch (const std::exception& e) {
            std::cerr << "Failed to load GA configuration: " << e.what() << std::endl;
            return ParameterSet();
        }
        
        std::random_device rd;
        std::mt19937 rng(rd());
        std::vector<ParameterSet> population = initPopulation(rng);
        
        ParameterSet bestSoFar;
        double bestFitnessSoFar = -1e9;

        int maxGenerations = config.getMaxGenerations();

        for (int gen = 0; gen < maxGenerations; gen++) {
            Racing::Manager rm(racingCfg); // Korzystanie z namespace
            auto results = rm.runRacing(population, [&](const ParameterSet &ps){
                return evaluator(ps);
            });

            for (auto &r : results) {
                if (r.fitness > bestFitnessSoFar) {
                    bestFitnessSoFar = r.fitness;
                    bestSoFar = r.parameterSet;
                }
            }

            population = breedNewPopulation(results);
        }
        
        return bestSoFar;
    }

private:
    MetaEAConfig config;
    Racing::Configuration racingCfg; // Korzystanie z namespace

    std::vector<ParameterSet> initPopulation(std::mt19937 &rng) {
        std::vector<ParameterSet> pop;
        pop.reserve(config.populationSize);

        std::uniform_int_distribution<int> distPop(50, 300);
        std::uniform_real_distribution<double> distMut(0.0, 1.0);

        for (int i = 0; i < config.populationSize; i++) {
            try {
                ParameterSet ps;
                ps.params["populationSize"] = std::to_string(distPop(rng));
                ps.params["mutationRate"] = std::to_string(distMut(rng));
                pop.push_back(ps);
            } catch (const std::exception& e) {
                std::cerr << "Error generating parameter set: " << e.what() << std::endl;
            }
        }
        return pop;
    }

    std::vector<ParameterSet> breedNewPopulation(const std::vector<TuningResult> &results)
    {
        auto sorted = results;
        std::sort(sorted.begin(), sorted.end(),
                  [](const TuningResult &a, const TuningResult &b){ return a.fitness > b.fitness; });
        
        std::vector<ParameterSet> newPop;
        newPop.reserve(config.populationSize);

        int eliteCount = std::max(2, config.populationSize / 5);
        for (int i = 0; i < eliteCount && i < static_cast<int>(sorted.size()); i++) {
            newPop.push_back(sorted[i].parameterSet);
        }
        
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dist(0, sorted.size()-1);

        while (static_cast<int>(newPop.size()) < config.populationSize && !sorted.empty()) {
            int idx1 = dist(rng);
            int idx2 = dist(rng);
            ParameterSet child = crossoverParams(sorted[idx1].parameterSet,
                                                 sorted[idx2].parameterSet);
            mutateParams(child, rng);
            newPop.push_back(child);
        }

        return newPop;
    }

    ParameterSet crossoverParams(const ParameterSet &p1, const ParameterSet &p2) {
        ParameterSet child;
        for (const auto &kv : p1.params) {
            if (p2.params.find(kv.first) != p2.params.end()) {
                child.params[kv.first] = (rand() % 2 == 0) ?
                    kv.second : p2.params.at(kv.first);
            } else {
                child.params[kv.first] = kv.second;
            }
        }
        for (const auto &kv : p2.params) {
            if (child.params.find(kv.first) == child.params.end()) {
                child.params[kv.first] = kv.second;
            }
        }
        return child;
    }

    void mutateParams(ParameterSet &ps, std::mt19937 &rng) {
        std::uniform_real_distribution<double> chanceDist(0.0, 1.0);
        if (ps.params.find("populationSize") != ps.params.end()) {
            double c = chanceDist(rng);
            if (c < 0.10) {
                int oldVal = std::stoi(ps.params["populationSize"]);
                double changeFactor = 0.9 + 0.2 * chanceDist(rng);
                int newVal = static_cast<int>(oldVal * changeFactor);
                if (newVal < 2) newVal = 2;
                ps.params["populationSize"] = std::to_string(newVal);
            }
        }

        if (ps.params.find("mutationRate") != ps.params.end()) {
            double c = chanceDist(rng);
            if (c < 0.10) {
                double oldVal = std::stod(ps.params["mutationRate"]);
                double delta = 0.05 * (chanceDist(rng) - 0.5) * 2.0;
                double newVal = oldVal + delta;
                if (newVal < 0.0) newVal = 0.0;
                if (newVal > 1.0) newVal = 1.0;
                ps.params["mutationRate"] = std::to_string(newVal);
            }
        }
    }
};

// [ADDED LINES FOR EXTENDED META-EA FROM THE ARTICLE]

// -------------------------------------------------
// Hybrid (1+lambda) ES + Racing demonstration
// -------------------------------------------------
class HybridOnePlusLambdaEA {
public:
    HybridOnePlusLambdaEA() : lambda(10) {}

    /**
     * Przykład metody "runHybridOnePlusLambdaEA" pokazującej ideę
     * łączenia (1+λ) ES z Racing:
     *  1) Mamy "aktualny" zestaw parametrów ps
     *  2) Generujemy λ potomków (wariacje parametrów)
     *  3) Oceniamy je przez Racing (minimalizując liczbę powtórzeń)
     *  4) Przyjmujemy najlepszego zamiast rodzica, jeśli lepszy
     */
    ParameterSet runHybridOnePlusLambdaEA(const ParameterSet &parent, const DNAInstance &exampleInstance)
    {
        // Implementacja pozostaje bez zmian, ale używamy Racing::Configuration
        // Bierzemy parent -> oceniamy (jeśli trzeba).
        TuningResult parentResult = evaluateParamSet(parent, exampleInstance);

        ParameterSet bestSoFar = parent;
        double bestFitSoFar = parentResult.fitness;

        // Definiujemy config do Racing
        Racing::Configuration rc; // Korzystanie z namespace
        rc.significanceLevel = 0.05;
        rc.maxTrialsPerCandidate = 10;
        rc.minTrialsBeforeElimination = 3;
        rc.useBootstrap = false;

        // Główna pętla – np. 5 iteracji
        for (int iteration = 0; iteration < 5; iteration++) {
            // Tworzymy populację kandydatów = parent + λ mutacji
            std::vector<ParameterSet> candidates;
            candidates.push_back(bestSoFar);
            for (int i = 0; i < lambda; i++) {
                candidates.push_back(mutateParameterSet(bestSoFar));
            }

            // Racing
            Racing::Manager rm(rc); // Korzystanie z namespace
            auto results = rm.runRacing(candidates, [&](const ParameterSet &ps){
                return evaluateParamSet(ps, exampleInstance);
            });

            // Szukamy najlepszego
            double localBestFit = -1e9;
            ParameterSet localBestSet;
            for (auto &r : results) {
                if (r.fitness > localBestFit) {
                    localBestFit = r.fitness;
                    localBestSet = r.parameterSet;
                }
            }

            // Jeśli lepszy od obecnego "bestSoFar" – aktualizujemy
            if (localBestFit > bestFitSoFar) {
                bestFitSoFar = localBestFit;
                bestSoFar = localBestSet;
            }
        }

        return bestSoFar;
    }

private:
    int lambda;

    // Przykładowa funkcja oceny "evaluateParamSet", uruchamia normalny GA
    // z danym zestawem parametrów i zwraca TuningResult.
    TuningResult evaluateParamSet(const ParameterSet &ps, const DNAInstance &exampleInstance)
    {
        // Get singleton instance
        auto& config = GAConfig::getInstance();
        
        try {
            config.loadFromFile("config.cfg");
        } catch (const std::exception& e) {
            std::cerr << "Failed to load GA configuration: " << e.what() << std::endl;
            TuningResult result;
            result.fitness = -std::numeric_limits<double>::infinity();
            return result;
        }
        
        // Update configuration using proper setters
        for (const auto& kv : ps.params) {
            if (kv.first == "populationSize") {
                config.setPopulationSize(std::stoi(kv.second));
            } else if (kv.first == "mutationRate") {
                config.setMutationRate(std::stod(kv.second));
            }
            // Note: crossoverType is handled through getCrossover below
        }
        
        // Create and run GA
        GeneticAlgorithm ga(
            config.getRepresentation(),
            config.getSelection(),
            config.getCrossover("order"), // Use default crossover type
            config.getMutation(),
            config.getReplacement(),
            config.getFitness(),
            config.getStopping(),
            config.getCache(),
            config
        );
        
        ga.run(exampleInstance);
        
        TuningResult result;
        result.parameterSet = ps;
        result.fitness = ga.getBestFitness();
        return result;
    }

    ParameterSet mutateParameterSet(const ParameterSet &parent)
    {
        ParameterSet child = parent;
        static std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> dist01(0.0, 1.0);

        if (child.params.find("populationSize") != child.params.end()) {
            double c = dist01(rng);
            if (c < 0.3) { // 30% szans
                int oldVal = std::stoi(child.params["populationSize"]);
                double fac = 0.8 + 0.4*dist01(rng); 
                int newVal = static_cast<int>(oldVal * fac);
                if (newVal < 2) newVal = 2;
                child.params["populationSize"] = std::to_string(newVal);
            }
        }
        if (child.params.find("mutationRate") != child.params.end()) {
            double c = dist01(rng);
            if (c < 0.3) {
                double oldVal = std::stod(child.params["mutationRate"]);
                double delta = 0.1*(dist01(rng) - 0.5)*2.0; 
                double newVal = oldVal + delta;
                if (newVal < 0.0) newVal = 0.0;
                if (newVal > 1.0) newVal = 1.0;
                child.params["mutationRate"] = std::to_string(newVal);
            }
        }
        return child;
    }
};

#endif //META_EA_H
