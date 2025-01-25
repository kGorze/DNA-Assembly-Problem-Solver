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
    // Initialize population
    auto population = m_representation->initializePopulation(m_config.populationSize);
    
    // Evaluate initial population
    evaluatePopulation(population, instance);
    updateGlobalBest(population, instance);
    logGenerationStats(population, instance, 0);

    // Main loop
    for (int generation = 1; generation <= m_config.maxGenerations; ++generation) {
        // Select parents for next generation
        std::vector<std::shared_ptr<Individual>> newPopulation;
        newPopulation.reserve(m_config.populationSize);

        while (newPopulation.size() < m_config.populationSize) {
            // Tournament selection
            std::uniform_int_distribution<> dis(0, population.size() - 1);
            int idx1 = dis(m_random->getGenerator());
            int idx2 = dis(m_random->getGenerator());
            
            auto parent1 = population[idx1]->getFitness() > population[idx2]->getFitness() ? 
                population[idx1] : population[idx2];
            
            idx1 = dis(m_random->getGenerator());
            idx2 = dis(m_random->getGenerator());
            
            auto parent2 = population[idx1]->getFitness() > population[idx2]->getFitness() ? 
                population[idx1] : population[idx2];

            // Crossover
            std::vector<int> child_genes;
            if (m_random->generateProbability() < m_config.crossoverProbability) {
                // Single-point crossover
                const auto& p1_genes = parent1->getGenes();
                const auto& p2_genes = parent2->getGenes();
                
                std::uniform_int_distribution<> crossover_point(0, p1_genes.size() - 1);
                int point = crossover_point(m_random->getGenerator());
                
                child_genes.insert(child_genes.end(), p1_genes.begin(), p1_genes.begin() + point);
                child_genes.insert(child_genes.end(), p2_genes.begin() + point, p2_genes.end());
            } else {
                child_genes = parent1->getGenes();
            }

            // Mutation
            if (m_random->generateProbability() < m_config.mutationProbability) {
                std::uniform_int_distribution<> pos(0, child_genes.size() - 1);
                int mutation_pos = pos(m_random->getGenerator());
                child_genes[mutation_pos] = 1 - child_genes[mutation_pos]; // Flip bit
            }

            auto child = std::make_shared<Individual>();
            child->setGenes(std::move(child_genes));
            newPopulation.push_back(child);
        }

        // Replace old population
        population = std::move(newPopulation);
        
        // Evaluate new population
        evaluatePopulation(population, instance);
        updateGlobalBest(population, instance);
        logGenerationStats(population, instance, generation);

        // Check termination condition
        if (m_globalBestFit >= m_config.targetFitness) {
            LOG_INFO("Target fitness reached in generation " + std::to_string(generation));
            break;
        }
    }

    return m_representation->toDNA(m_globalBestSolution);
}

void GeneticAlgorithm::evaluatePopulation(std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance) {
    for (auto& individual : population) {
        if (!individual) continue;
        
        if (!m_representation->isValid(individual, instance)) {
            individual->setFitness(-std::numeric_limits<double>::infinity());
            continue;
        }
        
        individual->setFitness(calculateFitness(individual, instance));
    }
}

void GeneticAlgorithm::updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance) {
    auto it = std::max_element(population.begin(), population.end(),
        [](const auto& a, const auto& b) { return a->getFitness() < b->getFitness(); });
    
    if (it != population.end() && (*it)->getFitness() > m_globalBestFit) {
        m_globalBestFit = (*it)->getFitness();
        m_globalBestSolution = *it;
    }
}

void GeneticAlgorithm::logGenerationStats(const std::vector<std::shared_ptr<Individual>>& population, const DNAInstance& instance, int generation) {
    double avgFitness = 0.0;
    double minFitness = std::numeric_limits<double>::infinity();
    double maxFitness = -std::numeric_limits<double>::infinity();
    
    for (const auto& individual : population) {
        if (!individual) continue;
        double fitness = individual->getFitness();
        if (std::isfinite(fitness)) {
            avgFitness += fitness;
            minFitness = std::min(minFitness, fitness);
            maxFitness = std::max(maxFitness, fitness);
        }
    }
    
    avgFitness /= population.size();
    
    std::stringstream ss;
    ss << "Generation " << generation 
       << " - Avg: " << avgFitness 
       << " Min: " << minFitness 
       << " Max: " << maxFitness 
       << " Global Best: " << m_globalBestFit;
    
    LOG_INFO(ss.str());
}

double GeneticAlgorithm::calculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) {
    if (!individual || !m_representation) return -std::numeric_limits<double>::infinity();
    
    // Convert individual's genes to DNA sequence
    std::string dna = m_representation->toDNA(individual);
    if (dna.empty()) return -std::numeric_limits<double>::infinity();
    
    // Calculate fitness based on:
    // 1. Length difference from target
    // 2. Number of matching k-mers with spectrum
    // 3. Connectivity between k-mers
    
    double lengthScore = 1.0 - std::abs(static_cast<double>(dna.length() - instance.getN())) / instance.getN();
    
    // Extract k-mers from DNA sequence
    std::vector<std::string> dnaKmers;
    for (size_t i = 0; i + instance.getK() <= dna.length(); ++i) {
        dnaKmers.push_back(dna.substr(i, instance.getK()));
    }
    
    // Count matching k-mers
    const auto& spectrum = instance.getSpectrum();
    int matches = 0;
    for (const auto& kmer : dnaKmers) {
        if (std::find(spectrum.begin(), spectrum.end(), kmer) != spectrum.end()) {
            matches++;
        }
    }
    double coverageScore = static_cast<double>(matches) / spectrum.size();
    
    // Calculate connectivity score
    double connectivityScore = 0.0;
    for (size_t i = 1; i < dnaKmers.size(); ++i) {
        const std::string& prev = dnaKmers[i-1];
        const std::string& curr = dnaKmers[i];
        if (prev.substr(1) == curr.substr(0, instance.getK()-1)) {
            connectivityScore += 1.0;
        }
    }
    connectivityScore /= std::max(1.0, static_cast<double>(dnaKmers.size() - 1));
    
    // Combine scores with weights
    double fitness = 0.4 * lengthScore + 0.3 * coverageScore + 0.3 * connectivityScore;
    return fitness;
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