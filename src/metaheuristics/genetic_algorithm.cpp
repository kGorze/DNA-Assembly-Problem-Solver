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

    std::mutex s_outputMutex;
}

GeneticAlgorithm::GeneticAlgorithm(
    std::unique_ptr<IRepresentation> representation,
    const GAConfig& config,
    bool debugMode)
    : m_representation(std::move(representation))
    , m_config(config)
    , m_debugMode(debugMode)
    , m_random(std::make_unique<Random>())
    , m_globalBestFit(-std::numeric_limits<double>::infinity())
{
    if (!m_representation) {
        throw std::invalid_argument("Invalid representation pointer");
    }
}

std::string GeneticAlgorithm::run(const DNAInstance& instance) {
    if (!m_representation) {
        LOG_ERROR("No representation set");
        return "";
    }

    try {
        // Initialize population
        initializePopulation(m_config.getPopulationSize(), instance);
        
        // Create stopping criteria
        auto stopping = m_config.getStopping();
        
        // Main loop
        int generation = 1;
        
        while (!stopping->shouldStop(generation, m_globalBestFit)) {
            // Selection
            std::vector<std::shared_ptr<Individual>> parents;
            parents.reserve(m_config.getPopulationSize());
            
            // Tournament selection
            for (int i = 0; i < m_config.getPopulationSize(); i++) {
                double bestFitness = -std::numeric_limits<double>::infinity();
                std::shared_ptr<Individual> selectedParent;
                
                for (int j = 0; j < m_config.getTournamentSize(); j++) {
                    int idx = m_random->getRandomInt(0, m_config.getPopulationSize() - 1);
                    double fitness = m_population[idx]->getFitness();
                    
                    if (fitness > bestFitness) {
                        bestFitness = fitness;
                        selectedParent = m_population[idx];
                    }
                }
                
                if (selectedParent) {
                    parents.push_back(selectedParent);
                }
            }
            
            // Crossover and Mutation
            std::vector<std::shared_ptr<Individual>> offspring;
            offspring.reserve(parents.size());
            
            for (size_t i = 0; i < parents.size() - 1; i += 2) {
                std::vector<std::shared_ptr<Individual>> parentPair = {parents[i], parents[i + 1]};
                
                if (m_random->generateProbability() < m_config.getCrossoverProbability()) {
                    auto crossover = m_config.getCrossover("");
                    auto children = crossover->crossover(parentPair, instance, std::shared_ptr<IRepresentation>(m_representation.get()));
                    offspring.insert(offspring.end(), children.begin(), children.end());
                } else {
                    offspring.push_back(parents[i]);
                    offspring.push_back(parents[i + 1]);
                }
            }
            
            // Handle odd number of parents
            if (parents.size() % 2 == 1) {
                offspring.push_back(parents.back());
            }
            
            // Mutation
            auto mutation = m_config.getMutation();
            for (auto& individual : offspring) {
                if (m_random->generateProbability() < m_config.getMutationRate()) {
                    mutation->mutate(individual, instance, std::shared_ptr<IRepresentation>(m_representation.get()));
                }
            }
            
            // Evaluate offspring
            evaluatePopulation(instance, offspring);
            
            // Replacement
            auto replacement = m_config.getReplacement();
            m_population = replacement->replace(m_population, offspring, instance, std::shared_ptr<IRepresentation>(m_representation.get()));
            
            // Update best solution
            updateGlobalBest(m_population, instance);
            
            // Log progress
            if (m_debugMode) {
                logGenerationStats(m_population, instance, generation);
            }
            
            generation++;
        }
        
        // Return best solution
        return m_bestDNA;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Error in genetic algorithm: " + std::string(e.what()));
        return "";
    }
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance& instance) {
    m_population.clear();
    m_population = m_representation->initializePopulation(popSize, instance);
    evaluatePopulation(instance, m_population);
    updateGlobalBest(m_population, instance);
}

void GeneticAlgorithm::evaluatePopulation(
    const DNAInstance& instance,
    const std::vector<std::shared_ptr<Individual>>& population)
{
    auto fitness = m_config.getFitness();
    for (auto& individual : population) {
        if (!individual) continue;
        double fitnessValue = fitness->calculateFitness(individual, instance, std::shared_ptr<IRepresentation>(m_representation.get()));
        individual->setFitness(fitnessValue);
    }
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<Individual>>& population,
    [[maybe_unused]] const DNAInstance& instance)
{
    for (const auto& individual : population) {
        if (!individual) continue;
        
        double fitness = individual->getFitness();
        if (fitness > m_globalBestFit) {
            m_globalBestFit = fitness;
            m_bestFitness = fitness;
            m_globalBestGenes = individual->getGenes();
            m_bestDNA = vectorToString(individual->getGenes());
        }
    }
}

void GeneticAlgorithm::logGenerationStats(
    const std::vector<std::shared_ptr<Individual>>& population,
    [[maybe_unused]] const DNAInstance& instance,
    int generation)
{
    double sumFitness = 0.0;
    double minFitness = std::numeric_limits<double>::infinity();
    double maxFitness = -std::numeric_limits<double>::infinity();
    int validCount = 0;
    
    for (const auto& individual : population) {
        if (!individual) continue;
        
        double fitness = individual->getFitness();
        sumFitness += fitness;
        minFitness = std::min(minFitness, fitness);
        maxFitness = std::max(maxFitness, fitness);
        validCount++;
    }
    
    double avgFitness = validCount > 0 ? sumFitness / validCount : 0.0;
    
    std::stringstream ss;
    ss << "Generation " << generation << ": ";
    ss << "Best=" << std::fixed << std::setprecision(4) << maxFitness << " ";
    ss << "Avg=" << avgFitness << " ";
    ss << "Min=" << minFitness;
    
    LOG_INFO(ss.str());
}

std::string GeneticAlgorithm::vectorToString(const std::vector<int>& vec) {
    if (vec.empty()) return "";
    
    std::stringstream ss;
    for (size_t i = 0; i < vec.size(); i++) {
        if (i > 0) ss << " ";
        ss << vec[i];
    }
    return ss.str();
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance& instance) {
    // Calculate theoretical maximum fitness based on instance parameters
    int k = instance.getK();
    int n = instance.getN();
    int spectrumSize = instance.getSpectrum().size();
    
    if (k <= 0 || n <= 0 || spectrumSize <= 0) {
        LOG_ERROR("Invalid instance parameters for theoretical fitness calculation");
        m_globalBestFit = 1.0;  // Default value
        return;
    }
    
    // Maximum possible score components
    double maxConnectivityScore = spectrumSize - 1;  // All k-mers perfectly connected
    double maxSpectrumCoverage = spectrumSize;       // All k-mers used
    double maxLengthScore = 1.0;                     // Perfect length match
    
    // Combine components with weights
    m_globalBestFit = maxConnectivityScore + maxSpectrumCoverage + maxLengthScore;
    
    LOG_INFO("Theoretical maximum fitness calculated: " + std::to_string(m_globalBestFit));
}

std::vector<std::vector<PreprocessedEdge>> GeneticAlgorithm::buildAdjacencyMatrix(
    const DNAInstance& instance) const {
    const auto& spectrum = instance.getSpectrum();
    const int k = instance.getK();
    const int deltaK = instance.getDeltaK();
    
    std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix(
        spectrum.size(),
        std::vector<PreprocessedEdge>(spectrum.size())
    );
    
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
                adjacencyMatrix[i][j] = PreprocessedEdge(j, weight, weight >= k - deltaK);
            }
        }
    }
    
    return adjacencyMatrix;
}

int GeneticAlgorithm::calculateEdgeWeight(const std::string& from, const std::string& to, int k) const {
    if (from.length() < static_cast<size_t>(k - 1) || to.length() < static_cast<size_t>(k)) return 0;
    
    std::string suffix = from.substr(from.length() - (k - 1));
    std::string prefix = to.substr(0, k - 1);
    
    return (suffix == prefix) ? 1 : 0;
}