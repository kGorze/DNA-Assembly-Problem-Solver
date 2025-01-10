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
                                         int generation) {
    if (generation % 10 != 0) return;

    // Calculate theoretical max fitness if not done yet
    if (m_theoreticalMaxFitness == 0.0) {
        calculateTheoreticalMaxFitness(instance);
    }

    double bestFit = -std::numeric_limits<double>::infinity();
    double worstFit = std::numeric_limits<double>::infinity();
    double sumFit = 0.0;

#pragma omp parallel for reduction(max:bestFit) reduction(min:worstFit) reduction(+:sumFit)
    for (size_t i = 0; i < pop.size(); i++) {
        // Use cache instead of direct calculation
        double fitVal = m_cache->getOrCalculateFitness(pop[i], instance, m_fitness, m_representation);
        bestFit = std::max(bestFit, fitVal);
        worstFit = std::min(worstFit, fitVal);
        sumFit += fitVal;
    }

    double avgFit = sumFit / pop.size();
    double relativeBestFit = (m_theoreticalMaxFitness > 0) ? 
        (bestFit / m_theoreticalMaxFitness) * 100.0 : 0.0;

    std::cout << "[Gen " << generation << "] "
              << "BestFit = " << bestFit 
              << " (" << std::fixed << std::setprecision(2) << relativeBestFit << "% of max "
              << m_theoreticalMaxFitness << ")"
              << ", WorstFit = " << worstFit
              << ", AvgFit = " << avgFit
              << ", PopSize = " << pop.size()
              << "\n";
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance &instance)
{
    population = m_representation->initializePopulation(popSize, instance);
}

void GeneticAlgorithm::updateGlobalBest(const std::vector<void*> &pop,
                                        const DNAInstance &instance)
{
    for (auto &ind : pop) {
        // Use cache instead of direct calculation
        double fitVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
        if (fitVal > m_globalBestFit) {
            m_globalBestFit = fitVal;

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
        
        if (generation % 1 == 0) {
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

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance &instance) {
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    if (spectrum.empty() || k <= 0) {
        m_theoreticalMaxFitness = 0.0;
        return;
    }

    // For OptimizedGraphBasedFitness, the theoretical maximum would be:
    // 1. Every edge in the path has weight 1 (best case)
    // 2. Every node is visited exactly once (perfect coverage)
    
    size_t n = spectrum.size();
    
    // Best possible edge score:
    // - In a path visiting all nodes once, we have (n-1) edges
    // - Best case: all edges have weight 1
    double bestEdgeScore = (n - 1);  // Each weight-1 edge contributes +1
    
    // Best possible coverage score:
    // - All nodes used exactly once (n unique nodes)
    // - No repeat usages
    double bestCoverageScore = n;  // n unique nodes, 0 repeats
    
    // Apply the weights from OptimizedGraphBasedFitness
    const double alpha = 0.7;  // Weight for edge score
    const double beta = 0.3;   // Weight for coverage score
    
    m_theoreticalMaxFitness = (alpha * bestEdgeScore) + (beta * bestCoverageScore);
    
    std::cout << "[MaxFitness] Theoretical maximum calculated:\n"
              << "- Best edge score: " << bestEdgeScore << " (weight: " << alpha << ")\n"
              << "- Best coverage score: " << bestCoverageScore << " (weight: " << beta << ")\n"
              << "- Total theoretical maximum: " << m_theoreticalMaxFitness << "\n";
}
