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
#include <iomanip>
#include <numeric>

std::mutex GeneticAlgorithm::outputMutex;


GeneticAlgorithm::GeneticAlgorithm(
    std::shared_ptr<IRepresentation> representation,
    std::shared_ptr<ISelection> selection,
    std::shared_ptr<ICrossover> crossover,
    std::shared_ptr<IMutation> mutation,
    std::shared_ptr<IReplacement> replacement,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IPopulationCache> cache,
    std::shared_ptr<IStopping> stopping,
    const GeneticConfig& config)
    : m_representation(std::move(representation))
    , m_selection(std::move(selection))
    , m_crossover(std::move(crossover))
    , m_mutation(std::move(mutation))
    , m_replacement(std::move(replacement))
    , m_fitness(std::move(fitness))
    , m_cache(std::move(cache))
    , m_stopping(std::move(stopping))
    , m_config(config)
{
    LOG_INFO("Initializing Genetic Algorithm with population size: " + std::to_string(m_config.getPopulationSize()));
}

GeneticAlgorithm::~GeneticAlgorithm() {
    population.clear();
    m_globalBestInd.reset();
}

void GeneticAlgorithm::logGenerationStats(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance,
    int generation)
{
    // Always calculate statistics
    if (m_theoreticalMaxFitness == 0.0) {
        calculateTheoreticalMaxFitness(instance);
    }

    double bestFit = -std::numeric_limits<double>::infinity();
    double avgFit = 0.0;
    double worstFit = std::numeric_limits<double>::infinity();
    int validSolutions = 0;
    
#pragma omp parallel for reduction(max:bestFit) reduction(+:avgFit) reduction(min:worstFit) reduction(+:validSolutions)
    for (size_t i = 0; i < pop.size(); i++) {
        if (!pop[i]) {
            LOG_ERROR("Null individual found at index " + std::to_string(i));
            continue;
        }
        
        double fitVal = pop[i]->getFitness();
        bestFit = std::max(bestFit, fitVal);
        avgFit += fitVal;
        worstFit = std::min(worstFit, fitVal);
        
        if (m_representation->isValid(pop[i], instance)) {
            validSolutions++;
        }
    }
    
    avgFit /= pop.size();
    
    double progress = bestFit / m_theoreticalMaxFitness * 100.0;
    
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << "Generation " << generation << ": ";
    ss << "Best=" << bestFit << " (" << progress << "%), ";
    ss << "Avg=" << avgFit << ", ";
    ss << "Worst=" << worstFit << ", ";
    ss << "Valid=" << validSolutions << "/" << pop.size();
    
    LOG_INFO(ss.str());
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance& instance)
{
    population = m_representation->initializePopulation(popSize, instance);
    
    int validCount = 0;
    int invalidCount = 0;
    int nullCount = 0;
    
    for (size_t i = 0; i < population.size(); i++) {
        if (!population[i]) {
            LOG_ERROR("Individual " + std::to_string(i) + " is nullptr");
            nullCount++;
            continue;
        }
        
        if(m_representation->isValid(population[i], instance)) {
            validCount++;
        } else {
            LOG_WARNING("Individual " + std::to_string(i) + " is invalid");
            invalidCount++;
            
            // Log invalid genes
            std::string invalidGenes;
            for(int gene : population[i]->getGenes()) {
                if(gene < 0 || gene >= (int)instance.getSpectrum().size()) {
                    invalidGenes += std::to_string(gene) + " ";
                }
            }
            if(!invalidGenes.empty()) {
                LOG_DEBUG("Invalid genes in individual " + std::to_string(i) + ": " + invalidGenes);
            }
        }
    }
    
    LOG_INFO("Population initialization complete:");
    LOG_INFO("  Valid individuals: " + std::to_string(validCount));
    LOG_INFO("  Invalid individuals: " + std::to_string(invalidCount));
    LOG_INFO("  Null individuals: " + std::to_string(nullCount));
    
    if(validCount == 0) {
        LOG_ERROR("No valid individuals in initial population!");
    }
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance)
{
    for (const auto& individual : pop) {
        if (!individual) continue;
        
        double fitness = individual->getFitness();
        
        if (fitness > m_globalBestFit) {
            m_globalBestFit = fitness;
            m_globalBestInd = std::make_shared<Individual>(*individual);
            m_config.setGlobalBestFitness(m_globalBestFit);
            
            // Convert best individual to DNA string
            auto dna = m_representation->toDNA(m_globalBestInd, instance);
            m_bestDNA = std::string(dna.begin(), dna.end());
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance& instance) {
    try {
        // Initialize random number generator
        std::random_device rd;
        generator = std::mt19937(rd());
        
        // Reset stopping criteria
        m_stopping->reset();
        
        // Clear population and cache
        population.clear();
        if (m_cache) {
            m_cache->clear();
        }
        
        // Initialize population
        initializePopulation(m_config.getPopulationSize(), instance);
        
        // Reserve cache space
        if (m_cache) {
            m_cache->reserve(m_config.getPopulationSize() * 2);  // Reserve space for parents and offspring
        }
        
        // Evaluate initial population
        for (auto& individual : population) {
            if (!individual) continue;
            double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
            individual->setFitness(fitness);
        }
        
        // Update best solution
        updateGlobalBest(population, instance);
        
        // Main loop
        int generation = 0;
        while (!m_stopping->shouldStop(generation, m_globalBestFit)) {
            // Log progress
            if (generation % m_config.getLogInterval() == 0) {
                logGenerationStats(population, instance, generation);
            }
            
            // Selection
            auto parents = m_selection->select(population, instance, m_fitness, m_representation);
            
            // Crossover
            auto offspring = m_crossover->crossover(parents, instance, m_representation);
            
            // Mutation
            for (auto& individual : offspring) {
                if (!individual) continue;
                m_mutation->mutate(individual, instance, m_representation);
            }
            
            // Remove invalid offspring
            offspring.erase(
                std::remove_if(offspring.begin(), offspring.end(),
                    [this, &instance](const auto& ind) {
                        return !ind || !m_representation->isValid(ind, instance);
                    }),
                offspring.end()
            );
            
            // Evaluate offspring
            for (auto& individual : offspring) {
                if (!individual) continue;
                double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
                individual->setFitness(fitness);
            }
            
            // Replacement
            population = m_replacement->replace(population, offspring, instance, m_representation);
            
            // Update best solution
            updateGlobalBest(population, instance);
            
            // Clean up cache if needed
            if (m_cache && generation % m_config.getCacheCleanupInterval() == 0) {
                m_cache->cleanup();
            }
            
            generation++;
        }
        
        // Final logging
        logGenerationStats(population, instance, generation);
        
        LOG_INFO("Genetic algorithm finished after " + std::to_string(generation) + " generations");
        LOG_INFO("Best fitness: " + std::to_string(m_globalBestFit));
        
        if (m_globalBestInd) {
            LOG_INFO("Best solution: " + m_bestDNA);
        } else {
            LOG_ERROR("No valid solution found!");
        }
        
    } catch (const std::exception& e) {
        LOG_ERROR("Genetic algorithm execution failed: " + std::string(e.what()));
        throw;
    }
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance& instance)
{
    // Calculate theoretical maximum fitness based on instance parameters
    int k = instance.getK();
    int n = instance.getN();
    int spectrumSize = instance.getSpectrum().size();
    
    if (k <= 0 || n <= 0 || spectrumSize <= 0) {
        LOG_ERROR("Invalid instance parameters for theoretical fitness calculation");
        m_theoreticalMaxFitness = 1.0;  // Default value
        return;
    }
    
    // Maximum possible score components
    double maxConnectivityScore = spectrumSize - 1;  // All k-mers perfectly connected
    double maxSpectrumCoverage = spectrumSize;       // All k-mers used
    double maxLengthScore = 1.0;                     // Perfect length match
    
    // Combine components with weights
    m_theoreticalMaxFitness = maxConnectivityScore + maxSpectrumCoverage + maxLengthScore;
    
    LOG_INFO("Theoretical maximum fitness calculated: " + std::to_string(m_theoreticalMaxFitness));
}

std::string GeneticAlgorithm::vectorToString(const std::vector<int>& vec)
{
    std::stringstream ss;
    for (int val : vec) {
        ss << val << " ";
    }
    return ss.str();
}

std::vector<std::vector<PreprocessedEdge>> GeneticAlgorithm::buildAdjacencyMatrix(
    const DNAInstance& instance) const
{
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    
    std::vector<std::vector<PreprocessedEdge>> matrix(
        spectrum.size(),
        std::vector<PreprocessedEdge>(spectrum.size())
    );
    
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
                if (weight > 0) {
                    matrix[i][j] = PreprocessedEdge(j, weight, true);
                }
            }
        }
    }
    
    return matrix;
}

int GeneticAlgorithm::calculateEdgeWeight(const std::string& from, const std::string& to, int k) const
{
    if (from.length() < k - 1 || to.length() < k - 1) return 0;
    
    // Compare suffix of 'from' with prefix of 'to'
    std::string suffix = from.substr(from.length() - (k - 1));
    std::string prefix = to.substr(0, k - 1);
    
    if (suffix == prefix) {
        return k - 1;  // Weight is the length of the overlap
    }
    
    return 0;
}

void GeneticAlgorithm::evolve(const DNAInstance& instance) {
    // Select parents
    auto parents = m_selection->select(population, instance, m_fitness, m_representation);
    
    if (parents.empty()) {
        LOG_ERROR("Selection returned empty parents list");
        return;
    }
    
    // Create offspring through crossover
    auto offspring = m_crossover->crossover(parents, instance, m_representation);
    
    if (offspring.empty()) {
        LOG_WARNING("Crossover produced no offspring, using parents instead");
        offspring = parents;
    }
    
    // Validate offspring before mutation
    offspring.erase(
        std::remove_if(offspring.begin(), offspring.end(),
            [](const auto& ind) { return !ind || ind->empty(); }),
        offspring.end()
    );
    
    if (offspring.empty()) {
        LOG_ERROR("No valid offspring after validation");
        return;
    }
    
    // Mutate offspring
    for (auto& child : offspring) {
        if (child && !child->empty()) {
            m_mutation->mutate(child, instance, m_representation);
        }
    }
    
    // Calculate fitness for population and offspring
    std::vector<double> populationFitness;
    populationFitness.reserve(population.size());
    for (const auto& individual : population) {
        if (individual && !individual->empty()) {
            populationFitness.push_back(m_cache->getOrCalculateFitness(individual, instance, m_fitness, m_representation));
        }
    }
    
    std::vector<double> offspringFitness;
    offspringFitness.reserve(offspring.size());
    for (const auto& individual : offspring) {
        if (individual && !individual->empty()) {
            offspringFitness.push_back(m_cache->getOrCalculateFitness(individual, instance, m_fitness, m_representation));
        }
    }
    
    // Replace old population with offspring
    population = m_replacement->replace(population, offspring, populationFitness, offspringFitness, instance, m_representation);
    
    // Validate final population
    if (population.empty()) {
        LOG_ERROR("Replacement produced empty population");
        population = parents;  // Fallback to parents if replacement fails
    }
}