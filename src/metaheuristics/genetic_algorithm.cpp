//
// Created by konrad_guest on 07/01/2025.
//

#include "metaheuristics/genetic_algorithm.h"
#include <iostream>
#include <random>
#include <chrono>

GeneticAlgorithm::GeneticAlgorithm(std::shared_ptr<IRepresentation> representation,
                                   std::shared_ptr<ISelection> selection,
                                   std::shared_ptr<ICrossover> crossover,
                                   std::shared_ptr<IMutation> mutation,
                                   std::shared_ptr<IReplacement> replacement,
                                   std::shared_ptr<IFitness> fitness,
                                   std::shared_ptr<IStopping> stopping)
    : m_representation(representation)
    , m_selection(selection)
    , m_crossover(crossover)
    , m_mutation(mutation)
    , m_replacement(replacement)
    , m_fitness(fitness)
    , m_stopping(stopping)
{}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance &instance)
{
    population = m_representation->initializePopulation(popSize, instance);
}

void GeneticAlgorithm::run(const DNAInstance &instance)
{
    // Example population size
    int POP_SIZE = 30; 
    initializePopulation(POP_SIZE, instance);

    int generation = 0;
    while (!m_stopping->stop(population, generation, instance, m_fitness, m_representation)) 
    {
        // 1) Selection
        auto parents = m_selection->select(population, instance, m_fitness, m_representation);
        // 2) Crossover
        auto offspring = m_crossover->crossover(parents, instance, m_representation);
        // 3) Mutation
        m_mutation->mutate(offspring, instance, m_representation);
        // 4) Replacement
        population = m_replacement->replace(population, offspring, instance, m_fitness, m_representation);

        generation++;
    }

    // Find best solution in final population
    double bestFit = -std::numeric_limits<double>::infinity();
    void* bestInd  = nullptr;

    for (auto &ind : population) {
        double fitVal = m_fitness->evaluate(ind, instance, m_representation);
        if (fitVal > bestFit) {
            bestFit = fitVal;
            bestInd = ind;
        }
    }

    // If you want to see what the best solution decodes to:
    if (bestInd) {
        // decode it via your representation
        std::string bestDNA = m_representation->decodeToDNA(bestInd, instance);
        std::cout << "[GA] End after generation " << generation 
                  << ", best fitness = " << bestFit
                  << ", length of best DNA = " << bestDNA.size() << "\n";
    }
    else {
        std::cout << "[GA] Population empty or no best individual found.\n";
    }
}