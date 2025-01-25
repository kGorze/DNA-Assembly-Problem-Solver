#include "../include/metaheuristics/genetic_algorithm.h"
#include "../include/configuration/genetic_algorithm_configuration.h"
#include "../include/metaheuristics/genetic_algorithm_runner.h"
#include "../include/utils/logging.h"
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <mutex>
#include <string>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <cmath>

namespace {
    // Helper function for safe string conversion
    template<typename T>
    std::string safeToString(const T& value) {
        try {
            return std::to_string(value);
        } catch (const std::exception&) {
            return "error";
        }
    }
}

std::mutex GeneticAlgorithm::s_outputMutex;

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
    , m_globalBestFit(-std::numeric_limits<double>::infinity())
    , m_globalBestInd(nullptr)
{
    if (!m_representation || !m_selection || !m_crossover || !m_mutation || 
        !m_replacement || !m_fitness || !m_cache || !m_stopping) {
        throw std::invalid_argument("All components must be non-null");
    }
    LOG_INFO("Initializing Genetic Algorithm with population size: " + safeToString(m_config.getPopulationSize()));
}

GeneticAlgorithm::~GeneticAlgorithm() {
    try {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_population.clear();
        m_globalBestInd = nullptr;
    } catch (const std::exception& e) {
        LOG_ERROR("Error during cleanup: " + std::string(e.what()));
    }
}

void GeneticAlgorithm::logGenerationStats(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance,
    int generation)
{
    if (pop.empty()) {
        LOG_ERROR("Empty population in generation " + safeToString(generation));
        return;
    }

    std::lock_guard<std::mutex> lock(s_outputMutex);

    // Always calculate statistics
    if (m_theoreticalMaxFitness == 0.0) {
        calculateTheoreticalMaxFitness(instance);
    }

    double bestFit = -std::numeric_limits<double>::infinity();
    double avgFit = 0.0;
    double worstFit = std::numeric_limits<double>::infinity();
    int validSolutions = 0;
    size_t nullCount = 0;
    
    for (const auto& individual : pop) {
        if (!individual) {
            nullCount++;
            continue;
        }
        
        const double fitVal = individual->getFitness();
        bestFit = std::max(bestFit, fitVal);
        avgFit += fitVal;
        worstFit = std::min(worstFit, fitVal);
        
        if (m_representation->isValid(individual, instance)) {
            validSolutions++;
        }
    }
    
    if (nullCount > 0) {
        LOG_WARNING("Found " + safeToString(nullCount) + " null individuals in generation " + safeToString(generation));
    }
    
    const size_t validPopSize = pop.size() - nullCount;
    if (validPopSize > 0) {
        avgFit /= validPopSize;
    }
    
    const double progress = (m_theoreticalMaxFitness > 0.0) ? 
        (bestFit / m_theoreticalMaxFitness * 100.0) : 0.0;
    
    std::stringstream ss;
    ss << std::fixed << std::setprecision(2);
    ss << "Generation " << generation << ": ";
    ss << "Best=" << bestFit << " (" << progress << "%), ";
    ss << "Avg=" << avgFit << ", ";
    ss << "Worst=" << worstFit << ", ";
    ss << "Valid=" << validSolutions << "/" << validPopSize;
    
    LOG_INFO(ss.str());

    // Call progress callback if set
    if (progressCallback) {
        progressCallback(generation, validSolutions, bestFit, avgFit, worstFit, progress);
    }
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance& instance)
{
    if (popSize <= 0) {
        throw std::invalid_argument("Population size must be positive");
    }

    std::lock_guard<std::mutex> lock(m_mutex);
    m_population = m_representation->initializePopulation(popSize, instance);
    
    if (m_population.empty()) {
        throw std::runtime_error("Failed to initialize population");
    }
    
    for (const auto& individual : m_population) {
        if (!individual || individual->empty()) {
            throw std::runtime_error("Invalid individual in initial population");
        }
    }
    
    int validCount = 0;
    int invalidCount = 0;
    int nullCount = 0;
    
    for (size_t i = 0; i < m_population.size(); i++) {
        if (!m_population[i]) {
            LOG_ERROR("Individual " + safeToString(i) + " is nullptr");
            nullCount++;
            continue;
        }
        
        if (m_representation->isValid(m_population[i], instance)) {
            validCount++;
        } else {
            LOG_WARNING("Individual " + safeToString(i) + " is invalid");
            invalidCount++;
            
            // Log invalid genes
            std::string invalidGenes;
            const auto& genes = m_population[i]->getGenes();
            for (const int gene : genes) {
                if (gene < 0 || gene >= static_cast<int>(instance.getSpectrum().size())) {
                    invalidGenes += safeToString(gene) + " ";
                }
            }
            if (!invalidGenes.empty()) {
                LOG_DEBUG("Invalid genes in individual " + safeToString(i) + ": " + invalidGenes);
            }
        }
    }
    
    LOG_INFO("Population initialization complete:");
    LOG_INFO("  Valid individuals: " + safeToString(validCount));
    LOG_INFO("  Invalid individuals: " + safeToString(invalidCount));
    LOG_INFO("  Null individuals: " + safeToString(nullCount));
    
    if (validCount == 0) {
        throw std::runtime_error("No valid individuals in initial population");
    }
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance)
{
    if (pop.empty()) {
        LOG_WARNING("Empty population in updateGlobalBest");
        return;
    }

    std::lock_guard<std::mutex> lock(m_mutex);
    for (const auto& individual : pop) {
        if (!individual) continue;
        
        try {
            double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
            if (fitness > m_globalBestFit) {
                m_globalBestFit = fitness;
                m_globalBestInd = std::make_shared<Individual>(*individual);
                m_config.setGlobalBestFitness(m_globalBestFit);
                
                // Convert best individual to DNA string
                const auto dna = m_representation->toDNA(m_globalBestInd, instance);
                m_bestDNA = std::string(dna.begin(), dna.end());
                
                LOG_INFO("New best solution found with fitness: " + safeToString(m_globalBestFit));
            }
        } catch (const std::exception& e) {
            LOG_WARNING("Error calculating fitness: " + std::string(e.what()));
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance& instance) {
    std::lock_guard<std::mutex> lock(m_mutex);
    
    try {
        // Reset stopping criteria
        m_stopping->reset();
        
        // Clear population and cache
        m_population.clear();
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
        for (auto& individual : m_population) {
            if (!individual) continue;
            try {
                const double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
                individual->setFitness(fitness);
            } catch (const std::exception& e) {
                LOG_ERROR("Failed to evaluate individual: " + std::string(e.what()));
            }
        }
        
        // Update best solution
        updateGlobalBest(m_population, instance);
        
        // Main loop
        int generation = 0;
        while (!m_stopping->shouldStop(generation, m_globalBestFit)) {
            try {
                // Log progress
                if (generation % m_config.getLogInterval() == 0) {
                    logGenerationStats(m_population, instance, generation);
                }
                
                // Selection
                auto parents = m_selection->select(m_population, instance, m_fitness, m_representation);
                
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
                    try {
                        const double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
                        individual->setFitness(fitness);
                    } catch (const std::exception& e) {
                        LOG_ERROR("Failed to evaluate offspring: " + std::string(e.what()));
                    }
                }
                
                // Replacement
                m_population = m_replacement->replace(m_population, offspring, instance, m_representation);
                
                // Update best solution
                updateGlobalBest(m_population, instance);
                
                // Clean up cache if needed
                if (m_cache && generation % m_config.getCacheCleanupInterval() == 0) {
                    m_cache->cleanup();
                }
                
                generation++;
            } catch (const std::exception& e) {
                LOG_ERROR("Error in generation " + safeToString(generation) + ": " + e.what());
                // Continue to next generation
            }
        }
        
        // Final logging
        logGenerationStats(m_population, instance, generation);
        
        LOG_INFO("Genetic algorithm finished after " + safeToString(generation) + " generations");
        LOG_INFO("Best fitness: " + safeToString(m_globalBestFit));
        
        if (m_globalBestInd) {
            LOG_INFO("Best solution: " + m_bestDNA);
        } else {
            throw std::runtime_error("No valid solution found");
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Fatal error in genetic algorithm: " + std::string(e.what()));
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