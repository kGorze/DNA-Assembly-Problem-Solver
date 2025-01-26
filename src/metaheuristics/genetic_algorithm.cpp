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

GeneticAlgorithm::GeneticAlgorithm(std::unique_ptr<IRepresentation> representation, const GAConfig& config)
    : m_representation(std::move(representation))
    , m_random(std::make_unique<Random>())
    , m_config(config)
    , m_globalBestFit(-std::numeric_limits<double>::infinity())
{
    if (!m_representation) {
        throw std::invalid_argument("Invalid representation");
    }
}

void GeneticAlgorithm::logGenerationStats(
    const std::vector<std::shared_ptr<Individual>>& pop,
    const DNAInstance& instance,
    int generation) {
    
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
            double fitness = calculateFitness(individual, instance);
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
}

void GeneticAlgorithm::initializePopulation(int popSize, const DNAInstance& instance) {
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

void GeneticAlgorithm::updateGlobalBest(const std::vector<std::shared_ptr<Individual>>& population, [[maybe_unused]] const DNAInstance& instance) {
    if (population.empty()) return;

    auto bestIndividual = *std::max_element(population.begin(), population.end(),
        [](const auto& a, const auto& b) { return a->getFitness() < b->getFitness(); });

    if (bestIndividual->getFitness() > m_globalBestFit) {
        m_globalBestFit = bestIndividual->getFitness();
        m_globalBestGenes = bestIndividual->getGenes();
        LOG_INFO("New best fitness: " + std::to_string(m_globalBestFit));
    }
}

std::string GeneticAlgorithm::run(const DNAInstance& instance) {
    try {
        // Initialize population
        m_population = m_representation->initializePopulation(m_config.getPopulationSize(), instance);
        if (m_population.empty()) {
            LOG_ERROR("Failed to initialize population");
            return "";
        }

        // Main loop
        for (int generation = 1; generation <= m_config.getMaxGenerations(); ++generation) {
            // Create new population
            std::vector<std::shared_ptr<Individual>> newPopulation;
            newPopulation.reserve(m_config.getPopulationSize());

            // Selection and crossover
            while (newPopulation.size() < static_cast<size_t>(m_config.getPopulationSize())) {
                // Tournament selection for parent 1
                std::shared_ptr<Individual> parent1;
                double bestFitness1 = -1e9;
                for (int j = 0; j < m_config.getTournamentSize(); ++j) {
                    int idx = m_random->getRandomInt(0, static_cast<int>(m_population.size() - 1));
                    double fitness = calculateFitness(m_population[idx], instance);
                    if (fitness > bestFitness1) {
                        bestFitness1 = fitness;
                        parent1 = m_population[idx];
                    }
                }

                // Tournament selection for parent 2
                std::shared_ptr<Individual> parent2;
                double bestFitness2 = -1e9;
                for (int j = 0; j < m_config.getTournamentSize(); ++j) {
                    int idx = m_random->getRandomInt(0, static_cast<int>(m_population.size() - 1));
                    double fitness = calculateFitness(m_population[idx], instance);
                    if (fitness > bestFitness2) {
                        bestFitness2 = fitness;
                        parent2 = m_population[idx];
                    }
                }

                // Crossover
                if (m_random->generateProbability() < m_config.getCrossoverProbability() &&
                    parent1 && parent2) {
                    auto child = std::make_shared<Individual>();
                    auto& parent1Genes = parent1->getGenes();
                    auto& parent2Genes = parent2->getGenes();

                    // Single-point crossover
                    int crossPoint = m_random->getRandomInt(0, static_cast<int>(parent1Genes.size() - 1));
                    std::vector<int> childGenes;
                    childGenes.insert(childGenes.end(), parent1Genes.begin(), parent1Genes.begin() + crossPoint);
                    childGenes.insert(childGenes.end(), parent2Genes.begin() + crossPoint, parent2Genes.end());
                    child->setGenes(childGenes);

                    // Mutation
                    if (m_random->generateProbability() < m_config.getMutationRate()) {
                        auto& genes = child->getGenes();
                        int pos1 = m_random->getRandomInt(0, static_cast<int>(genes.size() - 1));
                        int pos2 = m_random->getRandomInt(0, static_cast<int>(genes.size() - 1));
                        std::swap(genes[pos1], genes[pos2]);
                        child->setGenes(genes);
                    }

                    if (m_representation->isValid(child, instance)) {
                        newPopulation.push_back(child);
                    }
                }
            }

            // Update population
            m_population = std::move(newPopulation);

            // Update global best
            updateGlobalBest(m_population, instance);

            // Log progress
            logGenerationStats(m_population, instance, generation);

            // Check termination condition
            if (m_globalBestFit >= m_config.getTargetFitness()) {
                LOG_INFO("Target fitness reached");
                break;
            }
        }

        // Return best solution
        if (m_globalBestGenes.empty()) {
            LOG_ERROR("No valid solution found");
            return "";
        }

        return m_representation->toDNA(std::make_shared<Individual>(m_globalBestGenes), instance);

    } catch (const std::exception& e) {
        LOG_ERROR("Error in genetic algorithm: " + std::string(e.what()));
        return "";
    }
}

void GeneticAlgorithm::evaluatePopulation(const DNAInstance& instance, const std::vector<std::shared_ptr<Individual>>& population) {
    for (auto& individual : population) {
        if (!m_representation->isValid(individual, instance)) {
            individual->setFitness(-std::numeric_limits<double>::infinity());
            continue;
        }
        double fitness = calculateFitness(individual, instance);
        individual->setFitness(fitness);
    }
}

double GeneticAlgorithm::calculateFitness(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) {
    if (!individual || !m_representation) {
        return -std::numeric_limits<double>::infinity();
    }

    std::string dna = m_representation->toDNA(individual, instance);
    if (dna.empty()) {
        return -std::numeric_limits<double>::infinity();
    }

    // Length score - penalize differences from target length
    double lengthScore = 1.0 - std::abs(static_cast<double>(dna.length() - instance.getN())) / instance.getN();
    lengthScore = std::max(0.0, lengthScore);

    // Spectrum score - check how many k-mers from the spectrum are present
    double spectrumScore = 0.0;
    const auto& spectrum = instance.getSpectrum();
    if (!spectrum.empty()) {
        int matches = 0;
        for (size_t i = 0; i <= dna.length() - instance.getK(); ++i) {
            std::string kmer = dna.substr(i, instance.getK());
            if (std::find(spectrum.begin(), spectrum.end(), kmer) != spectrum.end()) {
                matches++;
            }
        }
        spectrumScore = static_cast<double>(matches) / spectrum.size();
    }

    // Connectivity score - check if nucleotides are properly connected
    double connectivityScore = 0.0;
    if (dna.length() >= 2) {
        int validConnections = 0;
        for (size_t i = 0; i < dna.length() - 1; ++i) {
            if (dna[i + 1] == 'A' || dna[i + 1] == 'C' || dna[i + 1] == 'G' || dna[i + 1] == 'T') {
                validConnections++;
            }
        }
        connectivityScore = static_cast<double>(validConnections) / (dna.length() - 1);
    }

    // Combine scores with weights
    const double lengthWeight = 0.4;
    const double spectrumWeight = 0.4;
    const double connectivityWeight = 0.2;

    return lengthWeight * lengthScore + 
           spectrumWeight * spectrumScore + 
           connectivityWeight * connectivityScore;
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

std::vector<std::vector<GeneticAlgorithm::PreprocessedEdge>> GeneticAlgorithm::buildAdjacencyMatrix(
    const DNAInstance& instance) const {
    const auto& spectrum = instance.getSpectrum();
    std::vector<std::vector<PreprocessedEdge>> matrix(spectrum.size());
    
    for (size_t i = 0; i < spectrum.size(); ++i) {
        for (size_t j = 0; j < spectrum.size(); ++j) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], instance.getK());
                if (weight > 0) {
                    matrix[i].push_back({static_cast<int>(j), weight});
                }
            }
        }
    }
    
    return matrix;
}

int GeneticAlgorithm::calculateEdgeWeight(const std::string& from, const std::string& to, int k) const {
    if (from.length() < static_cast<size_t>(k - 1) || to.length() < static_cast<size_t>(k)) return 0;
    
    std::string suffix = from.substr(from.length() - (k - 1));
    std::string prefix = to.substr(0, k - 1);
    
    return (suffix == prefix) ? 1 : 0;
}