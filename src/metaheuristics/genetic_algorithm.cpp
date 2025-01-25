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

GeneticAlgorithm::GeneticAlgorithm(const GeneticConfig& config, std::unique_ptr<IRepresentation> representation)
    : m_config(config)
    , m_representation(std::move(representation))
    , m_random(std::make_unique<Random>()) {}

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

void GeneticAlgorithm::updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance) {
    if (population.empty()) return;

    auto bestIndividual = *std::min_element(population.begin(), population.end(),
        [](const auto& a, const auto& b) { return a->getFitness() < b->getFitness(); });

    if (!m_globalBestInd || bestIndividual->getFitness() < m_globalBestFit) {
        m_globalBestFit = bestIndividual->getFitness();
        m_globalBestInd = std::make_shared<Individual>(*bestIndividual);
        m_bestDNA = m_representation->toString(bestIndividual, instance);
        LOG_INFO("New best fitness: " + std::to_string(m_globalBestFit));
        LOG_INFO("Best DNA: " + m_bestDNA);
    }
}

std::string GeneticAlgorithm::run(const DNAInstance& instance) {
    if (!m_representation) {
        throw std::runtime_error("Representation not set");
    }

    // Initialize population
    m_population = m_representation->initializePopulation(m_config.populationSize);
    
    // Evaluate initial population
    evaluatePopulation();
    
    // Main loop
    for (size_t generation = 0; generation < m_config.maxGenerations; ++generation) {
        // Selection
        auto parents = selectParents();
        
        // Crossover
        auto offspring = crossover(parents);
        
        // Mutation
        for (auto& individual : offspring) {
            if (m_random->generateProbability() < m_config.mutationProbability) {
                mutate(individual);
            }
        }
        
        // Replace population
        m_population = std::move(offspring);
        
        // Evaluate new population
        evaluatePopulation();
        
        // Update best solution
        updateBestSolution();
        
        // Log progress
        if (generation % 100 == 0) {
            LOG_INFO("Generation " << generation << ": Best fitness = " << m_bestFitness);
        }
        
        // Check termination condition
        if (m_bestFitness >= m_config.targetFitness) {
            LOG_INFO("Target fitness reached at generation " << generation);
            break;
        }
    }
    
    return m_representation->toDNA(m_bestSolution);
}

void GeneticAlgorithm::evaluatePopulation() {
    for (auto& individual : m_population) {
        if (!m_representation->isValid(individual)) {
            individual->setFitness(0.0);
            individual->setValid(false);
            continue;
        }
        
        double fitness = calculateFitness(individual);
        individual->setFitness(fitness);
        individual->setValid(true);
    }
}

std::vector<std::shared_ptr<Individual>> GeneticAlgorithm::selectParents() {
    std::vector<std::shared_ptr<Individual>> parents;
    parents.reserve(m_config.populationSize);
    
    // Tournament selection
    std::uniform_int_distribution<size_t> dist(0, m_population.size() - 1);
    
    while (parents.size() < m_config.populationSize) {
        size_t idx1 = dist(m_random->getGenerator());
        size_t idx2 = dist(m_random->getGenerator());
        
        if (m_population[idx1]->getFitness() > m_population[idx2]->getFitness()) {
            parents.push_back(m_population[idx1]);
        } else {
            parents.push_back(m_population[idx2]);
        }
    }
    
    return parents;
}

std::vector<std::shared_ptr<Individual>> GeneticAlgorithm::crossover(
    const std::vector<std::shared_ptr<Individual>>& parents) {
    std::vector<std::shared_ptr<Individual>> offspring;
    offspring.reserve(parents.size());
    
    for (size_t i = 0; i < parents.size(); i += 2) {
        if (i + 1 >= parents.size()) {
            offspring.push_back(std::make_shared<Individual>(*parents[i]));
            continue;
        }
        
        if (m_random->generateProbability() < m_config.crossoverProbability) {
            auto [child1, child2] = performCrossover(parents[i], parents[i + 1]);
            offspring.push_back(std::move(child1));
            offspring.push_back(std::move(child2));
        } else {
            offspring.push_back(std::make_shared<Individual>(*parents[i]));
            offspring.push_back(std::make_shared<Individual>(*parents[i + 1]));
        }
    }
    
    return offspring;
}

std::pair<std::shared_ptr<Individual>, std::shared_ptr<Individual>> GeneticAlgorithm::performCrossover(
    const std::shared_ptr<Individual>& parent1,
    const std::shared_ptr<Individual>& parent2) {
    auto child1 = std::make_shared<Individual>();
    auto child2 = std::make_shared<Individual>();
    
    const auto& genes1 = parent1->getGenes();
    const auto& genes2 = parent2->getGenes();
    
    std::uniform_int_distribution<size_t> dist(0, genes1.size() - 1);
    size_t crossoverPoint = dist(m_random->getGenerator());
    
    std::vector<int> childGenes1;
    std::vector<int> childGenes2;
    childGenes1.reserve(genes1.size());
    childGenes2.reserve(genes2.size());
    
    // First child gets first part from parent1, second part from parent2
    childGenes1.insert(childGenes1.end(), genes1.begin(), genes1.begin() + crossoverPoint);
    childGenes1.insert(childGenes1.end(), genes2.begin() + crossoverPoint, genes2.end());
    
    // Second child gets first part from parent2, second part from parent1
    childGenes2.insert(childGenes2.end(), genes2.begin(), genes2.begin() + crossoverPoint);
    childGenes2.insert(childGenes2.end(), genes1.begin() + crossoverPoint, genes1.end());
    
    child1->setGenes(std::move(childGenes1));
    child2->setGenes(std::move(childGenes2));
    
    return {child1, child2};
}

void GeneticAlgorithm::mutate(std::shared_ptr<Individual>& individual) {
    auto genes = individual->getGenes();
    std::uniform_int_distribution<size_t> posDist(0, genes.size() - 1);
    std::uniform_int_distribution<int> valDist(0, 3);
    
    size_t pos = posDist(m_random->getGenerator());
    genes[pos] = valDist(m_random->getGenerator());
    
    individual->setGenes(std::move(genes));
}

void GeneticAlgorithm::updateBestSolution() {
    auto it = std::max_element(m_population.begin(), m_population.end(),
        [](const auto& a, const auto& b) {
            return a->getFitness() < b->getFitness();
        });
    
    if (it != m_population.end() && (*it)->getFitness() > m_bestFitness) {
        m_bestSolution = *it;
        m_bestFitness = (*it)->getFitness();
    }
}

double GeneticAlgorithm::calculateFitness(const std::shared_ptr<Individual>& individual) {
    // TODO: Implement fitness calculation based on the problem requirements
    return 0.0;
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

std::string GeneticAlgorithm::vectorToString(const std::vector<int>& vec)
{
    std::stringstream ss;
    for (int val : vec) {
        ss << val << " ";
    }
    return ss.str();
}

std::vector<std::vector<PreprocessedEdge>> GeneticAlgorithm::buildAdjacencyMatrix(
    const DNAInstance& instance) const {
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

int GeneticAlgorithm::calculateEdgeWeight(const std::string& from, const std::string& to, int k) const {
    if (from.empty() || to.empty()) return 0;
    
    auto fromLen = static_cast<int>(from.length());
    auto toLen = static_cast<int>(to.length());
    
    if (fromLen < k - 1 || toLen < k - 1) return 0;

    std::string fromSuffix = from.substr(fromLen - k + 1);
    std::string toPrefix = to.substr(0, k - 1);
    
    return fromSuffix == toPrefix ? k - 1 : 0;
}