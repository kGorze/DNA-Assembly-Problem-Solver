//
// Created by konrad_guest on 23/01/2025.
//

#ifndef META_EA_H
#define META_EA_H

#include "racing.h"  // Dodaj ten nagłówek
#include "../generator/dna_generator.h"
#include "configuration/genetic_algorithm_configuration.h" // Dołączenie GAConfig
#include "metaheuristics/genetic_algorithm.h"  
#include "metaheuristics/representation.h"
#include "dna/dna_instance.h"
#include "parameters.h"
#include <memory>
#include <vector>
#include <random>
#include <chrono>

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

class IMetaEA {
public:
    virtual ~IMetaEA() = default;
    virtual TuningResult evaluateParamSet(const ParameterSet& params, const DNAInstance& instance) = 0;
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
        // Create local config
        GAConfig config;
        if (!config.loadFromFile("config.cfg")) {
            throw std::runtime_error("Failed to load configuration");
        }
        
        std::random_device rd;
        std::mt19937 rng(rd());
        std::vector<ParameterSet> population = initPopulation(rng);
        
        ParameterSet bestSoFar;
        double bestFitnessSoFar = -1e9;

        for (int gen = 0; gen < config.getPopulationSize(); gen++) {  // Use population size as a proxy for generations
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

    TuningResult evaluateParamSet(const ParameterSet& params, const DNAInstance& instance) {
        // Create local config
        GAConfig config;
        
        // Set parameters from ParameterSet using setters
        config.setPopulationSize(params.getInt("populationSize"));
        config.setMutationRate(params.getDouble("mutationRate"));
        config.setCrossoverProbability(params.getDouble("crossoverRate"));
        config.setTournamentSize(params.getInt("tournamentSize"));
        
        // Create result object
        TuningResult result;
        result.parameterSet = params;
        result.fitness = 0.0;
        result.executionTime = 0.0;
        
        try {
            auto start = std::chrono::high_resolution_clock::now();
            
            // Create representation
            auto representation = std::make_unique<PermutationRepresentation>();

            // Create genetic algorithm
            GeneticAlgorithm ga(std::move(representation), config);

            // Run the algorithm
            auto solution = ga.run(instance);
            
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = end - start;
            
            result.fitness = calculateFitness(solution, instance);
            result.executionTime = diff.count();
            result.extraMetrics["converged"] = 1.0; // Success
        } catch (const std::exception& e) {
            result.extraMetrics["converged"] = 0.0; // Failure
            LOG_ERROR("Error in evaluateParamSet: " + std::string(e.what()));
        }
        
        return result;
    }

private:
    MetaEAConfig config;
    Racing::Configuration racingCfg; // Korzystanie z namespace

    double calculateFitness(const std::string& solution, const DNAInstance& instance) const {
        // Calculate Hamming distance between solution and target DNA
        const std::string& targetDNA = instance.getDNA();
        if (solution.length() != targetDNA.length()) {
            return 0.0;  // Invalid solution
        }
        
        int matches = 0;
        for (size_t i = 0; i < solution.length(); ++i) {
            if (solution[i] == targetDNA[i]) {
                matches++;
            }
        }
        
        return static_cast<double>(matches) / solution.length();
    }

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
class HybridOnePlusLambdaEA : public IMetaEA {
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
        // Evaluate parent first
        TuningResult parentResult = evaluateParamSet(parent, exampleInstance);
        
        ParameterSet bestSoFar = parent;
        double bestFitSoFar = parentResult.fitness;

        // Define racing config
        Racing::Configuration rc;
        rc.significanceLevel = 0.05;
        rc.maxTrialsPerCandidate = 10;
        rc.minTrialsBeforeElimination = 3;
        rc.useBootstrap = false;

        // Main loop - 5 iterations
        for (int iteration = 0; iteration < 5; iteration++) {
            // Create candidate population = parent + λ mutations
            std::vector<ParameterSet> candidates;
            candidates.push_back(bestSoFar);
            for (int i = 0; i < lambda; i++) {
                candidates.push_back(mutateParameterSet(bestSoFar));
            }

            // Run racing
            Racing::Manager rm(rc);
            auto results = rm.runRacing(candidates, [&](const ParameterSet &ps) {
                return evaluateParamSet(ps, exampleInstance);
            });

            // Update best if found better
            for (const auto &result : results) {
                if (result.fitness > bestFitSoFar) {
                    bestFitSoFar = result.fitness;
                    bestSoFar = result.parameterSet;
                }
            }
        }
        
        return bestSoFar;
    }

    TuningResult evaluateParamSet(const ParameterSet& params, const DNAInstance& instance) override {
        try {
            // Create local config
            GAConfig config;
            
            // Set parameters from ParameterSet using setters
            config.setPopulationSize(params.getInt("populationSize"));
            config.setMutationRate(params.getDouble("mutationRate"));
            config.setCrossoverProbability(params.getDouble("crossoverRate"));
            config.setTournamentSize(params.getInt("tournamentSize"));
            
            // Create representation
            auto representation = std::make_unique<PermutationRepresentation>();

            // Create genetic algorithm
            GeneticAlgorithm ga(std::move(representation), config);

            // Run the algorithm
            auto solution = ga.run(instance);
            
            // Return result
            TuningResult result;
            result.parameterSet = params;
            result.fitness = calculateFitness(solution, instance);
            result.executionTime = 0.0;  // TODO: Add timing
            return result;
            
        } catch (const std::exception& e) {
            LOG_ERROR("Error in HybridOnePlusLambdaEA::evaluateParamSet: " + std::string(e.what()));
            return TuningResult{params, 0.0, 0.0};
        }
    }

private:
    int lambda;

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

    double calculateFitness(const std::string& solution, const DNAInstance& instance) {
        return instance.calculateFitness(solution);
    }
};

#endif //META_EA_H
