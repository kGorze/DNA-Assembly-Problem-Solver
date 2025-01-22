//
// Created by konrad_guest on 07/01/2025.
// SMART


#include "metaheuristics/genetic_algorithm.h"
#include <iostream>
#include <random>
#include <chrono>

std::mutex GeneticAlgorithm::outputMutex;

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
{ }

GeneticAlgorithm::~GeneticAlgorithm() {
    // Shared pointers automatycznie zwolnią pamięć,
    // więc nie ma potrzeby wywoływania `delete`.
    population.clear();
    m_globalBestInd.reset();
}

void GeneticAlgorithm::logGenerationStats(
    const std::vector<std::shared_ptr<std::vector<int>>>& pop,
    const DNAInstance& instance,
    int generation)
{
    if (generation % 10 != 0) return;

    if (m_theoreticalMaxFitness == 0.0) {
        calculateTheoreticalMaxFitness(instance);
    }

    double bestFit = -std::numeric_limits<double>::infinity();
    double avgFit = 0.0;
    
#pragma omp parallel for reduction(max:bestFit) reduction(+:avgFit)
    for (size_t i = 0; i < pop.size(); i++) {
        double fitVal = m_cache->getOrCalculateFitness(pop[i], instance, m_fitness, m_representation);
        if (fitVal > bestFit) bestFit = fitVal;
        avgFit += fitVal;
    }
    
    avgFit /= pop.size();
    double progress = (generation * 100.0) / GAConfig::getInstance().getMaxGenerations();
    double relativeBestFit = (m_theoreticalMaxFitness > 0) ? 
        (bestFit / m_theoreticalMaxFitness) * 100.0 : 0.0;

    // Format status with underscores instead of colons
    std::ostringstream status;
    status << "Best" << "_" << std::fixed << std::setprecision(2) << relativeBestFit 
           << "% Avg" << "_" << avgFit;

    // Synchronized output
    std::lock_guard<std::mutex> lock(outputMutex);
    std::cout << "PROGRESS_UPDATE:" << m_processId << ":" 
              << progress << ":" << status.str() << ":" 
              << relativeBestFit << std::endl;
    std::cout.flush();
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance &instance)
{
    population = m_representation->initializePopulation(popSize, instance);
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<std::vector<int>>> &pop,
    const DNAInstance &instance)
{
    for (const auto &ind : pop) {
        if (!ind) continue;
        double fitVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
        if (fitVal > m_globalBestFit) {
            m_globalBestFit = fitVal;
            m_globalBestInd = std::make_shared<std::vector<int>>(*ind);
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance &instance)
{
    PROFILE_FUNCTION();
    population = m_representation->initializePopulation(GAConfig::getInstance().getPopulationSize(), instance);
    m_cache->updatePopulation(population, instance, m_fitness, m_representation);

    int generation = 0;
    int maxGenerations = GAConfig::getInstance().getMaxGenerations();
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    std::vector<std::shared_ptr<std::vector<int>>> parents;
    
    while (!m_stopping->stop(population, generation, instance, m_fitness, m_representation)) {
        
        // 1) Selection
        {
            PROFILE_SCOPE("selection");
            parents = m_selection->select(population, instance, m_fitness, m_representation);
        }
        
        // 2) Crossover
        {
            PROFILE_SCOPE("crossover"); 
            offspring = m_crossover->crossover(parents, instance, m_representation);
        }
        
        // 3) Mutation
        {
            PROFILE_SCOPE("mutation");
            m_mutation->mutate(offspring, instance, m_representation);
        }
        
        // 4) Cache update
        {
            PROFILE_SCOPE("cache_update");
            m_cache->updatePopulation(offspring, instance, m_fitness, m_representation);
        }
        
        // 5) Replacement
        {
            PROFILE_SCOPE("replacement");
            population = m_replacement->replace(population, offspring, instance, m_fitness, m_representation);
        }

        // Update best fitness
        double currentBestFitness = -std::numeric_limits<double>::infinity();
        {
            for (auto& ind : population) {
                double fitnessVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
                if (fitnessVal > currentBestFitness) {
                    currentBestFitness = fitnessVal;
                }
            }
            
            auto adaptiveCrossover = std::dynamic_pointer_cast<AdaptiveCrossover>(m_crossover);
            if (adaptiveCrossover) {
                adaptiveCrossover->updateFeedback(currentBestFitness);
            }
        }
        
        if (generation % 5 == 0) {
            PROFILE_SCOPE("best");
            updateGlobalBest(population, instance);
        }
        
        if (generation % 1 == 0) {
            PROFILE_SCOPE("stats");
            logGenerationStats(population, instance, generation);
        }
        
        // Call progress callback with current status
        if (progressCallback) {
            progressCallback(generation, maxGenerations, currentBestFitness);
        }
        
        generation++;
    }

    // Final update
    updateGlobalBest(population, instance);
    if (m_globalBestInd) {
        m_bestDNA = m_representation->decodeToDNA(m_globalBestInd, instance);
        
        // Send final progress update
        if (progressCallback) {
            progressCallback(generation, maxGenerations, m_globalBestFit);
        }
        
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

    size_t n = spectrum.size();
    double bestEdgeScore = (n - 1);  
    double bestCoverageScore = n;  
    
    const double alpha = 0.7;  
    const double beta = 0.3;   
    
    m_theoreticalMaxFitness = (alpha * bestEdgeScore) + (beta * bestCoverageScore);
    
    std::cout << "[MaxFitness] Theoretical maximum calculated:\n"
              << "- Best edge score: " << bestEdgeScore << " (weight: " << alpha << ")\n"
              << "- Best coverage score: " << bestCoverageScore << " (weight: " << beta << ")\n"
              << "- Total theoretical maximum: " << m_theoreticalMaxFitness << "\n";
}
