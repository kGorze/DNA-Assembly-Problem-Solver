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
        
        double fitVal = m_cache->getOrCalculateFitness(pop[i], instance, m_fitness, m_representation);
        if (fitVal > bestFit) bestFit = fitVal;
        if (fitVal < worstFit) worstFit = fitVal;
        avgFit += fitVal;
        
        if (m_representation->isValid(pop[i], instance)) {
            validSolutions++;
        }
    }
    
    avgFit /= pop.size();
    double progress = (generation * 100.0) / m_config.getMaxGenerations();
    double relativeBestFit = (m_theoreticalMaxFitness > 0) 
        ? (bestFit / m_theoreticalMaxFitness) * 100.0
        : 0.0;

    // Always output detailed statistics in debug mode
    LOG_INFO("Generation " + std::to_string(generation) + " statistics:");
    LOG_INFO("  Progress: " + std::to_string(progress) + "%");
    LOG_INFO("  Best fitness: " + std::to_string(bestFit) + " (" + std::to_string(relativeBestFit) + "%)");
    LOG_INFO("  Average fitness: " + std::to_string(avgFit));
    LOG_INFO("  Worst fitness: " + std::to_string(worstFit));
    LOG_INFO("  Valid solutions: " + std::to_string(validSolutions) + "/" + std::to_string(pop.size()));
    
    if (m_globalBestInd && m_globalBestFit > -std::numeric_limits<double>::infinity()) {
        LOG_INFO("  Global best fitness: " + std::to_string(m_globalBestFit));
        
        // Analyze best path if available
        auto pathAnalyzer = std::make_shared<PathAnalyzer>();
        auto bestPath = *m_globalBestInd;
        auto adjacencyMatrix = buildAdjacencyMatrix(instance);
        auto analysis = pathAnalyzer->analyzePath(bestPath, adjacencyMatrix);
        
        LOG_INFO("  Best path analysis:");
        LOG_INFO("    Unique nodes used: " + std::to_string(analysis.uniqueNodesUsed));
        LOG_INFO("    Edges weight 1: " + std::to_string(analysis.edgesWeight1));
        LOG_INFO("    Edges weight 2/3: " + std::to_string(analysis.edgesWeight2or3));
        
        // Additional path analysis
        std::stringstream pathStr;
        pathStr << "  Path: [";
        for (size_t i = 0; i < bestPath.size(); ++i) {
            if (i > 0) pathStr << ", ";
            pathStr << bestPath[i];
        }
        pathStr << "]";
        LOG_DEBUG(pathStr.str());
        
        // Log spectrum fragments used
        LOG_DEBUG("  Spectrum fragments used:");
        const auto& spectrum = instance.getSpectrum();
        for (int idx : bestPath) {
            if (idx >= 0 && idx < static_cast<int>(spectrum.size())) {
                LOG_DEBUG("    " + std::to_string(idx) + ": " + spectrum[idx]);
            }
        }
    }

    // Standard progress update
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
    LOG_DEBUG("Initializing population with size " + std::to_string(popSize));
    
    population = m_representation->initializePopulation(popSize, instance);
    
    // Validate population
    int validCount = 0;
    int nullCount = 0;
    int invalidCount = 0;
    
    for(size_t i = 0; i < population.size(); i++) {
        if(!population[i]) {
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
            for(auto gene : *population[i]) {
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

void GeneticAlgorithm::run(const DNAInstance& instance) 
{
    LOG_INFO("Starting genetic algorithm run");
    LOG_DEBUG("Instance parameters:");
    LOG_DEBUG("  Spectrum size: " + std::to_string(instance.getSpectrum().size()));
    LOG_DEBUG("  k: " + std::to_string(instance.getK()));
    LOG_DEBUG("  deltaK: " + std::to_string(instance.getDeltaK()));
    LOG_DEBUG("  N: " + std::to_string(instance.getN()));
    LOG_DEBUG("  L_Neg: " + std::to_string(instance.getLNeg()));
    LOG_DEBUG("  L_Poz: " + std::to_string(instance.getLPoz()));
    
    // Validate configuration parameters
    if (m_config.getPopulationSize() <= 0) {
        LOG_ERROR("Invalid population size: " + std::to_string(m_config.getPopulationSize()));
        throw std::runtime_error("Population size must be positive");
    }
    
    if (m_config.getMaxGenerations() <= 0) {
        LOG_ERROR("Invalid max generations: " + std::to_string(m_config.getMaxGenerations()));
        throw std::runtime_error("Max generations must be positive");
    }
    
    if (m_config.getMutationRate() < 0.0 || m_config.getMutationRate() > 1.0) {
        LOG_ERROR("Invalid mutation rate: " + std::to_string(m_config.getMutationRate()));
        throw std::runtime_error("Mutation rate must be between 0 and 1");
    }
    
    if (m_config.getCrossoverProbability() < 0.0 || m_config.getCrossoverProbability() > 1.0) {
        LOG_ERROR("Invalid crossover probability: " + std::to_string(m_config.getCrossoverProbability()));
        throw std::runtime_error("Crossover probability must be between 0 and 1");
    }
    
    try {
        // Initialize population
        LOG_INFO("Initializing population with size " + std::to_string(m_config.getPopulationSize()));
        initializePopulation(m_config.getPopulationSize(), instance);
        
        // Validate initial population
        if (population.empty()) {
            LOG_ERROR("Failed to initialize population - population is empty");
            throw std::runtime_error("Population initialization failed");
        }
        
        int validIndividuals = 0;
        for (const auto& ind : population) {
            if (ind && m_representation->isValid(ind, instance)) {
                validIndividuals++;
            }
        }
        
        if (validIndividuals == 0) {
            LOG_ERROR("No valid individuals in initial population");
            throw std::runtime_error("Population initialization failed - no valid individuals");
        }
        
        LOG_INFO("Initial population created with " + std::to_string(validIndividuals) + " valid individuals");
        
        // Initial evaluation
        LOG_INFO("Performing initial population evaluation");
        updateGlobalBest(population, instance);
        LOG_INFO("Initial best fitness: " + std::to_string(m_globalBestFit));
        
        int generation = 0;
        int stagnationCounter = 0;
        double previousBestFitness = m_globalBestFit;
        bool shouldStop = false;
        
        while (!shouldStop && generation < m_config.getMaxGenerations()) {
            try {
                LOG_DEBUG("Starting generation " + std::to_string(generation));
                
                // Selection
                LOG_DEBUG("Performing selection");
                auto parents = m_selection->select(population, instance, m_fitness, m_representation);
                if (parents.empty()) {
                    LOG_ERROR("Selection returned empty parents list in generation " + std::to_string(generation));
                    throw std::runtime_error("Selection failed - no parents selected");
                }
                LOG_DEBUG("Selected " + std::to_string(parents.size()) + " parents");
                
                // Crossover
                LOG_DEBUG("Performing crossover");
                auto offspring = m_crossover->crossover(parents, instance, m_representation);
                if (offspring.empty()) {
                    LOG_WARNING("Crossover produced no offspring, using parents as offspring");
                    offspring = parents;
                }
                LOG_DEBUG("Working with " + std::to_string(offspring.size()) + " offspring");
                
                // Mutation
                LOG_DEBUG("Performing mutation");
                int mutatedCount = 0;
                for (auto& child : offspring) {
                    if (!child || child->empty()) {
                        continue;
                    }
                    try {
                        if ((double)rand() / RAND_MAX < m_config.getMutationRate()) {
                            m_mutation->mutate(child, instance, m_representation);
                            mutatedCount++;
                        }
                    } catch (const std::exception& e) {
                        LOG_ERROR("Mutation failed for offspring: " + std::string(e.what()));
                    }
                }
                LOG_DEBUG("Mutated " + std::to_string(mutatedCount) + " offspring");
                
                // Calculate fitness
                std::vector<double> parentFitness;
                std::vector<double> offspringFitness;
                
                LOG_DEBUG("Calculating fitness for parents and offspring");
                for (const auto& ind : population) {
                    try {
                        parentFitness.push_back(m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation));
                    } catch (const std::exception& e) {
                        LOG_ERROR("Failed to calculate parent fitness: " + std::string(e.what()));
                        parentFitness.push_back(-std::numeric_limits<double>::infinity());
                    }
                }
                
                for (const auto& ind : offspring) {
                    try {
                        offspringFitness.push_back(m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation));
                    } catch (const std::exception& e) {
                        LOG_ERROR("Failed to calculate offspring fitness: " + std::string(e.what()));
                        offspringFitness.push_back(-std::numeric_limits<double>::infinity());
                    }
                }
                
                // Replacement
                LOG_DEBUG("Performing replacement");
                try {
                    auto newPopulation = m_replacement->replace(population, offspring, parentFitness, offspringFitness, instance, m_representation);
                    if (newPopulation.empty()) {
                        LOG_WARNING("Replacement returned empty population, keeping current population");
                    } else {
                        population = std::move(newPopulation);
                    }
                    LOG_DEBUG("Population size after replacement: " + std::to_string(population.size()));
                } catch (const std::exception& e) {
                    LOG_ERROR("Replacement failed: " + std::string(e.what()));
                    LOG_WARNING("Keeping previous population due to replacement failure");
                }
                
                // Update statistics and log
                updateGlobalBest(population, instance);
                logGenerationStats(population, instance, generation);
                
                // Check for improvement
                if (m_globalBestFit <= previousBestFitness) {
                    stagnationCounter++;
                    if (stagnationCounter >= 50) { // 50 generations without improvement
                        LOG_INFO("Algorithm stagnated for 50 generations");
                        shouldStop = true;
                    }
                } else {
                    stagnationCounter = 0;
                    previousBestFitness = m_globalBestFit;
                }
                
                // Check stopping condition
                if (!shouldStop) {
                    shouldStop = m_stopping->stop(population, instance, generation, m_globalBestFit);
                }
                
                if (shouldStop) {
                    LOG_INFO("Stopping condition met at generation " + std::to_string(generation));
                    LOG_INFO("Final best fitness: " + std::to_string(m_globalBestFit));
                }
                
            } catch (const std::exception& e) {
                LOG_ERROR("Error in generation " + std::to_string(generation) + ": " + std::string(e.what()));
                if (++generation >= m_config.getMaxGenerations()) {
                    LOG_INFO("Maximum generations reached after error");
                    break;
                }
                continue;
            }
            
            generation++;
        }
        
        LOG_INFO("Genetic algorithm completed");
        LOG_INFO("Total generations: " + std::to_string(generation));
        LOG_INFO("Best fitness achieved: " + std::to_string(m_globalBestFit));
        
        if (m_globalBestInd) {
            LOG_DEBUG("Best solution path length: " + std::to_string(m_globalBestInd->size()));
            LOG_DEBUG("Best solution: " + vectorToString(*m_globalBestInd));
            
            // Validate final solution
            if (!m_representation->isValid(m_globalBestInd, instance)) {
                LOG_WARNING("Best solution is not valid according to representation!");
            }
        } else {
            LOG_WARNING("No valid solution found during the run");
        }
        
    } catch (const std::exception& e) {
        LOG_ERROR("Fatal error in genetic algorithm: " + std::string(e.what()));
        throw;
    }
}

std::string GeneticAlgorithm::vectorToString(const std::vector<int>& vec) {
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) ss << ", ";
        ss << vec[i];
    }
    ss << "]";
    return ss.str();
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance& instance) {
    LOG_DEBUG("Calculating theoretical maximum fitness");
    
    int n = instance.getN();
    int k = instance.getK();
    int pathLength = n - k + 1;
    
    // Maximum possible values for each component
    double maxConnectivityScore = 1.0;  // All edges weight 1
    double maxUniquenessScore = 1.0;    // All nodes used exactly once
    double maxRepetitionScore = 1.0;    // No repeated nodes
    
    // Weights from fitness calculation
    constexpr double CONNECTIVITY_WEIGHT = 0.5;
    constexpr double UNIQUENESS_WEIGHT = 0.3;
    constexpr double REPETITION_WEIGHT = 0.2;
    
    m_theoreticalMaxFitness = (CONNECTIVITY_WEIGHT * maxConnectivityScore) +
                             (UNIQUENESS_WEIGHT * maxUniquenessScore) +
                             (REPETITION_WEIGHT * maxRepetitionScore);
    
    LOG_DEBUG("Theoretical maximum fitness calculated: " + std::to_string(m_theoreticalMaxFitness));
}

int GeneticAlgorithm::calculateEdgeWeight(const std::string& from, const std::string& to, int k) const {
    if (from.length() < k - 1 || to.length() < k - 1) {
        return 0;
    }
    
    // Check overlap between suffix of 'from' and prefix of 'to'
    std::string fromSuffix = from.substr(from.length() - (k - 1));
    std::string toPrefix = to.substr(0, k - 1);
    
    if (fromSuffix == toPrefix) {
        return 1;  // Perfect overlap
    }
    
    // Check for partial overlaps
    int maxOverlap = 0;
    for (int i = 1; i < k - 1; ++i) {
        if (fromSuffix.substr(i) == toPrefix.substr(0, k - 1 - i)) {
            maxOverlap = k - 1 - i;
        }
    }
    
    return maxOverlap > 0 ? 2 : 0;  // Weight 2 for partial overlaps
}

std::vector<std::vector<PreprocessedEdge>> GeneticAlgorithm::buildAdjacencyMatrix(const DNAInstance& instance) const {
    const auto& spectrum = instance.getSpectrum();
    int size = spectrum.size();
    std::vector<std::vector<PreprocessedEdge>> matrix(size, std::vector<PreprocessedEdge>(size));
    
    LOG_DEBUG("Building adjacency matrix for spectrum size: " + std::to_string(size));
    
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (i != j) {
                const std::string& from = spectrum[i];
                const std::string& to = spectrum[j];
                int weight = calculateEdgeWeight(from, to, instance.getK());
                matrix[i][j] = PreprocessedEdge{weight > 0, weight};
            }
        }
    }
    
    return matrix;
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