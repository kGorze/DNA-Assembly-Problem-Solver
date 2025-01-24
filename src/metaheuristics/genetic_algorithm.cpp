#include "../include/metaheuristics/genetic_algorithm.h"
#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../include/metaheuristics/genetic_algorithm_runner.h"
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <mutex>
#include <string>
#include "../include/utils/logging.h"
#include <sstream>

std::mutex GeneticAlgorithm::outputMutex;


GeneticAlgorithm::GeneticAlgorithm(
    std::shared_ptr<IRepresentation> representation,
    std::shared_ptr<ISelection> selection,
    std::shared_ptr<ICrossover> crossover,
    std::shared_ptr<IMutation> mutation,
    std::shared_ptr<IReplacement> replacement,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IStopping> stopping,
    std::shared_ptr<IPopulationCache> cache,
    GAConfig& config
) : m_representation(representation),
    m_selection(selection),
    m_crossover(crossover),
    m_mutation(mutation),
    m_replacement(replacement),
    m_fitness(fitness),
    m_stopping(stopping),
    m_cache(cache),
    m_config(config),
    m_globalBestFit(-std::numeric_limits<double>::infinity())
{
    LOG_INFO("Initializing Genetic Algorithm with population size: " + std::to_string(m_config.getPopulationSize()));
}

GeneticAlgorithm::~GeneticAlgorithm() {
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
    double progress = (generation * 100.0) / m_config.getMaxGenerations();
    double relativeBestFit = (m_theoreticalMaxFitness > 0) 
        ? (bestFit / m_theoreticalMaxFitness) * 100.0
        : 0.0;

    std::ostringstream status;
    status << "Best_" << std::fixed << std::setprecision(2) << relativeBestFit
           << "% Avg_" << avgFit;

    std::lock_guard<std::mutex> lock(outputMutex);
    std::cout << "PROGRESS_UPDATE:" << m_processId << ":"
              << progress << ":"
              << status.str() << ":"
              << relativeBestFit
              << std::endl;
    std::cout.flush();
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance &instance)
{
    population = m_representation->initializePopulation(popSize, instance);

    // **Debug Log**: Sprawdzenie populacji
    std::cout << "[DEBUG GA] Population initialized with size: " << population.size() << "\n";
    for(size_t i = 0; i < population.size(); i++) {
        auto &ind = population[i];
        if(!ind){
            std::cerr << "[ERROR GA] Individual " << i << " is nullptr\n";
            continue;
        }
        for(auto gene : *ind){
            if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                std::cerr << "[ERROR GA] Individual " << i << " has out-of-range gene: " << gene << "\n";
            }
        }
    }
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<std::vector<int>>>& pop,
    const DNAInstance& instance)
{
    for (const auto& individual : pop) {
        if (!individual) continue;
        
        double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
        
        if (fitness > m_globalBestFit) {
            m_globalBestFit = fitness;
            m_globalBestInd = std::make_shared<std::vector<int>>(*individual);
            m_config.setGlobalBestFitness(m_globalBestFit);
            
            // Convert best individual to DNA string
            std::stringstream ss;
            for (int gene : *m_globalBestInd) {
                ss << gene;
            }
            m_bestDNA = ss.str();
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance& instance) {
    // Initialize population
    population = m_representation->initializePopulation(m_config.getPopulationSize(), instance);
    
    if (population.empty()) {
        LOG_ERROR("Failed to initialize population");
        return;
    }
    
    int generation = 0;
    updateGlobalBest(population, instance);
    
    // Main evolution loop
    while (!m_stopping->stop(population, instance, generation, m_globalBestFit)) {
        evolve(instance);
        
        // Update global best
        updateGlobalBest(population, instance);
        
        // Log progress
        logGenerationStats(population, instance, generation);
        
        generation++;
    }
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance &instance)
{
    // Obliczamy teoretyczne maksymalne fitness
    // Dla prostego fitness to po prostu długość DNA
    // Dla bardziej złożonych implementacji można nadpisać tę metodę
    m_theoreticalMaxFitness = instance.getSize();
}

void GeneticAlgorithm::evolve(const DNAInstance& instance) {
    // Select parents
    auto parents = m_selection->select(population, instance, m_fitness, m_representation);
    
    // Create offspring through crossover
    auto offspring = m_crossover->crossover(parents, instance, m_representation);
    
    // Mutate offspring
    for (auto& child : offspring) {
        m_mutation->mutate(child, instance, m_representation);
    }
    
    // Calculate fitness for population and offspring
    std::vector<double> populationFitness;
    populationFitness.reserve(population.size());
    for (const auto& individual : population) {
        populationFitness.push_back(m_cache->getOrCalculateFitness(individual, instance, m_fitness, m_representation));
    }
    
    std::vector<double> offspringFitness;
    offspringFitness.reserve(offspring.size());
    for (const auto& individual : offspring) {
        offspringFitness.push_back(m_cache->getOrCalculateFitness(individual, instance, m_fitness, m_representation));
    }
    
    // Replace old population with offspring
    population = m_replacement->replace(population, offspring, populationFitness, offspringFitness, instance, m_representation);
}

double GeneticAlgorithm::getBestFitness() const {
    return m_globalBestFit;
}

std::shared_ptr<std::vector<int>> GeneticAlgorithm::getBestIndividual() const {
    return m_globalBestInd;
}