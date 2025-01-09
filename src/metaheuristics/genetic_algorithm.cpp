//
// Created by konrad_guest on 07/01/2025.
//

#include "metaheuristics/genetic_algorithm.h"
#include <iostream>
#include <random>
#include <chrono>

GeneticAlgorithm::GeneticAlgorithm(
    std::shared_ptr<IRepresentation> representation,
    std::shared_ptr<ISelection> selection,
    std::shared_ptr<ICrossover> crossover,
    std::shared_ptr<IMutation> mutation,
    std::shared_ptr<IReplacement> replacement,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IStopping> stopping,
    std::shared_ptr<IPopulationCache> cache)
    : m_representation(representation)
    , m_selection(selection)
    , m_crossover(crossover)
    , m_mutation(mutation)
    , m_replacement(replacement)
    , m_fitness(fitness)
    , m_stopping(stopping)
    , m_cache(cache)
    , m_globalBestInd(nullptr)
    , m_globalBestFit(-std::numeric_limits<double>::infinity())
{}

void GeneticAlgorithm::logGenerationStats(const std::vector<void*> &pop,
                                          const DNAInstance &instance,
                                          int generation)
{
    // Inicjuj zmienne do przechowywania statystyk
    double bestFit    = -std::numeric_limits<double>::infinity();
    double worstFit   =  std::numeric_limits<double>::infinity();
    double sumFit     = 0.0;
    int count         = 0;
    void* bestInd     = nullptr;

    // Iteracja po populacji
    for (auto &ind : pop) {
        double fitVal = m_fitness->evaluate(ind, instance, m_representation);
        sumFit  += fitVal;
        count++;

        if (fitVal > bestFit) {
            bestFit = fitVal;
            bestInd = ind;
        }
        if (fitVal < worstFit) {
            worstFit = fitVal;
        }
    }

    double avgFit = (count > 0) ? (sumFit / count) : 0.0;

    // Wypisz statystyki
    std::cout << "[Gen " << generation << "] "
              << "BestFit = "   << bestFit 
              << ", WorstFit = " << worstFit
              << ", AvgFit = "   << avgFit
              << ", PopSize = "  << pop.size()
              << "\n";
    
    // if (bestInd) {
    //     std::string bestDNA = m_representation->decodeToDNA(bestInd, instance);
    //     std::cout << "   (Best DNA length: " << bestDNA.size() << ")\n";
    // }
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance &instance)
{
    population = m_representation->initializePopulation(popSize, instance);
}

void GeneticAlgorithm::updateGlobalBest(const std::vector<void*> &pop,
                                        const DNAInstance &instance)
{
    for (auto &ind : pop) {
        double fitVal = m_fitness->evaluate(ind, instance, m_representation);
        if (fitVal > m_globalBestFit) {
            m_globalBestFit = fitVal;

            // SKOPIUJ osobnika, by uniknąć problemu z usuwaniem oryginalnego w replacement
            if(m_globalBestInd) {
                delete static_cast<std::vector<int>*>(m_globalBestInd);
            }
            auto orig = static_cast<std::vector<int>*>(ind);
            m_globalBestInd = new std::vector<int>(*orig);
        }
    }
}


void GeneticAlgorithm::run(const DNAInstance &instance)
{
    PROFILE_FUNCTION();
    population = m_representation->initializePopulation(GAConfig::getInstance().getPopulationSize(), instance);
    m_cache->updatePopulation(population, instance, m_fitness, m_representation);

    int generation = 0;
    std::vector<void*> offspring;
    std::vector<void*> parents;
    
    while (!m_stopping->stop(population, generation, instance, m_fitness, m_representation)) {
        {
            PROFILE_SCOPE("selection");
            parents = m_selection->select(population, instance, m_fitness, m_representation);
        }
        
        {
            PROFILE_SCOPE("crossover"); 
            offspring = m_crossover->crossover(parents, instance, m_representation);
        }
        
        {
            PROFILE_SCOPE("mutation");
            m_mutation->mutate(offspring, instance, m_representation);
        }
        
        {
            PROFILE_SCOPE("cache_update");
            m_cache->updatePopulation(offspring, instance, m_fitness, m_representation);
        }
        
        {
            PROFILE_SCOPE("replacement");
            population = m_replacement->replace(population, offspring, instance, m_fitness, m_representation);
        }
        
        if (generation % 5 == 0) {
            PROFILE_SCOPE("best");
            updateGlobalBest(population, instance);
        }
        
        if (generation % 25 == 0) {
            PROFILE_SCOPE("stats");
            logGenerationStats(population, instance, generation);
        }
        
        generation++;
    }

    // Po zakończeniu pętli - używamy m_globalBestInd
    if (m_globalBestInd) {
        m_bestDNA = m_representation->decodeToDNA(m_globalBestInd, instance);
        std::cout << "[GA] End after generation " << generation 
                  << ", best-ever fitness = " << m_globalBestFit
                  << ", length of best DNA = " << m_bestDNA.size() << "\n"
                  << "Cache size: " << m_cache->size() << "\n";
    }
}