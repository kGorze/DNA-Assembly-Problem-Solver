#include "../../include/metaheuristics/genetic_algorithm.h"
#include "../../include/configuration/genetic_algorithm_configuration.h"
#include "../../include/metaheuristics/genetic_algorithm_runner.h"
#include "../../include/utils/logging.h"
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

// Initialize static mutex
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
    std::lock_guard<std::mutex> lock(m_mutex);
    m_population.clear();
    m_globalBestInd = nullptr;
}

void GeneticAlgorithm::logGenerationStats(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance,
    int generation) {
    
    std::lock_guard<std::mutex> lock(m_mutex);
    
    if (pop.empty()) {
        LOG_WARNING("Empty population in generation " + std::to_string(generation));
        return;
    }

    double sumFitness = 0.0;
    double minFitness = std::numeric_limits<double>::infinity();
    double maxFitness = -std::numeric_limits<double>::infinity();
    int validCount = 0;

    for (const auto& individual : pop) {
        if (!individual) continue;
        
        try {
            double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
            sumFitness += fitness;
            minFitness = std::min(minFitness, fitness);
            maxFitness = std::max(maxFitness, fitness);
            if (individual->isValid()) validCount++;
        } catch (const std::exception& e) {
            LOG_WARNING("Error calculating fitness: " + std::string(e.what()));
        }
    }

    double avgFitness = sumFitness / pop.size();
    double validRatio = static_cast<double>(validCount) / pop.size();

    LOG_INFO("Generation " + std::to_string(generation) + 
             " stats: min=" + std::to_string(minFitness) + 
             ", avg=" + std::to_string(avgFitness) + 
             ", max=" + std::to_string(maxFitness) + 
             ", valid=" + std::to_string(validRatio));

    if (progressCallback) {
        progressCallback(generation, 
                       static_cast<int>(pop.size()),
                       avgFitness,
                       maxFitness,
                       m_globalBestFit,
                       m_theoreticalMaxFitness);
    }
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance& instance) {
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
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance) {
    
    std::lock_guard<std::mutex> lock(m_mutex);
    
    for (const auto& individual : pop) {
        if (!individual) continue;
        
        try {
            double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
            if (fitness > m_globalBestFit) {
                m_globalBestFit = fitness;
                m_globalBestInd = std::make_shared<Individual>(*individual);
                m_bestDNA = m_representation->toString(individual->getGenes());
            }
        } catch (const std::exception& e) {
            LOG_WARNING("Error calculating fitness: " + std::string(e.what()));
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance& instance) {
    std::lock_guard<std::mutex> lock(m_mutex);
    
    try {
        // Clear any existing state
        m_population.clear();
        m_globalBestFit = -std::numeric_limits<double>::infinity();
        m_globalBestInd = nullptr;
        
        // Initialize population
        initializePopulation(m_config.getPopulationSize(), instance);
        
        // Reserve space in cache
        m_cache->reserve(m_config.getPopulationSize() * 2);
        
        int generation = 0;
        while (!m_stopping->shouldStop(generation, m_globalBestFit)) {
            // Evaluate fitness for all individuals
            for (auto& individual : m_population) {
                if (!individual) continue;
                
                try {
                    double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
                    individual->setFitness(fitness);
                } catch (const std::exception& e) {
                    LOG_WARNING("Error calculating fitness: " + std::string(e.what()));
                    continue;
                }
            }
            
            // Update global best
            updateGlobalBest(m_population, instance);
            
            // Log progress
            logGenerationStats(m_population, instance, generation);
            
            // Selection
            auto parents = m_selection->select(m_population, m_config.getParentCount());
            
            // Crossover
            auto offspring = m_crossover->crossover(parents, instance, m_representation);
            
            // Mutation
            for (auto& individual : offspring) {
                if (!individual) continue;
                m_mutation->mutate(individual, instance, m_representation);
            }
            
            // Evaluate offspring fitness
            for (auto& individual : offspring) {
                if (!individual) continue;
                
                try {
                    double fitness = m_fitness->calculateFitness(individual, instance, m_representation);
                    individual->setFitness(fitness);
                } catch (const std::exception& e) {
                    LOG_WARNING("Error calculating offspring fitness: " + std::string(e.what()));
                    continue;
                }
            }
            
            // Replacement
            m_population = m_replacement->replace(m_population, offspring);
            
            // Clear cache periodically
            if (generation % 10 == 0) {
                m_cache->clear();
            }
            
            ++generation;
        }
        
        LOG_INFO("Genetic algorithm completed after " + std::to_string(generation) + " generations");
        LOG_INFO("Best fitness: " + std::to_string(m_globalBestFit));
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in genetic algorithm: " + std::string(e.what()));
        throw;
    }
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance& instance)
{
    std::lock_guard<std::mutex> lock(m_mutex);
    
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
    
    // Pre-allocate the matrix with the correct size
    std::vector<std::vector<PreprocessedEdge>> matrix;
    matrix.reserve(spectrum.size());
    for (size_t i = 0; i < spectrum.size(); ++i) {
        matrix.emplace_back(spectrum.size());
    }
    
    // Build the matrix
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

int GeneticAlgorithm::calculateEdgeWeight(
    const std::string& from,
    const std::string& to,
    int k) const
{
    if (from.length() < k - 1 || to.length() < k - 1) return 0;
    
    // Compare suffix of 'from' with prefix of 'to'
    std::string_view suffix(from.data() + from.length() - (k - 1), k - 1);
    std::string_view prefix(to.data(), k - 1);
    
    return suffix == prefix ? k - 1 : 0;
}
}