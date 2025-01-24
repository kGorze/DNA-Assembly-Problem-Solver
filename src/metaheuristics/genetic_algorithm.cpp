#include "metaheuristics/genetic_algorithm.h"
#include "configuration/genetic_algorithm_configuration.h" // żeby mieć dostęp do GAConfig
#include "metaheuristics/genetic_algorithm_runner.h"  // for updateConfigWithInstanceParams
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <mutex>
#include <string>
#include "utils/logging.h"

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
    LOG_INFO("Initializing Genetic Algorithm with population size: " + std::to_string(m_config.populationSize));
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
    const std::vector<std::shared_ptr<std::vector<int>>> &pop,
    const DNAInstance &instance)
{
    for (const auto &ind : pop) {
        if (!ind) {
            std::cerr << "[ERROR GA] Best candidate is nullptr\n";
            continue;
        }
        double fitVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
        if (fitVal > m_globalBestFit) {
            m_globalBestFit = fitVal;
            m_globalBestInd = std::make_shared<std::vector<int>>(*ind);
            // **Debug Log**: Nowy globalny best
            std::cout << "[DEBUG GA] New global best fitness: " << m_globalBestFit << "\n";
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance &instance)
{
    LOG_INFO("Starting Genetic Algorithm optimization");
    DEBUG_LOG("Instance size: " + std::to_string(instance.getSize()));
    
    auto startTime = std::chrono::high_resolution_clock::now();
    
    // Log the maxGenerations value at the start
    std::cout << "[DEBUG GA] Starting with maxGenerations = " << m_config.getMaxGenerations() << std::endl;
    
    population = m_representation->initializePopulation(m_config.populationSize, instance);
    m_cache->updatePopulation(population, instance, m_fitness, m_representation);

    // **Debug Log**: Po inicjalizacji populacji
    std::cout << "[DEBUG GA] After initialization, population size: " << population.size() << "\n";
    for(size_t i = 0; i < population.size(); i++) {
        auto &ind = population[i];
        if(!ind){
            std::cerr << "[ERROR GA] Individual " << i << " is nullptr after initialization\n";
            continue;
        }
        for(auto gene : *ind){
            if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                std::cerr << "[ERROR GA] Individual " << i << " has out-of-range gene after initialization: " << gene << "\n";
            }
        }
    }

    int generation = 0;
    // Always get maxGenerations from config when needed instead of storing locally
    std::cout << "[DEBUG GA] Current maxGenerations = " << m_config.getMaxGenerations() << std::endl;

    // Z góry obliczamy teoretyczne maksymalne fitness
    calculateTheoreticalMaxFitness(instance);

    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    std::vector<std::shared_ptr<std::vector<int>>> parents;
    
    std::cout << "[DEBUG GA] Starting main GA loop with maxGenerations = " << m_config.getMaxGenerations() << "\n";

    while (!m_stopping->stop(population, generation, instance, m_fitness, m_representation)) {
        std::cout << "[DEBUG GA] Generation " << generation << " of " << m_config.getMaxGenerations() << "\n";
        
        {
            PROFILE_SCOPE("selection");
            parents = m_selection->select(population, instance, m_fitness, m_representation);
        }

        // **Debug Log**: Po selekcji
        std::cout << "[DEBUG GA] After selection, parents size: " << parents.size() << "\n";
        for(size_t i = 0; i < parents.size(); i++) {
            auto &ind = parents[i];
            if(!ind){
                std::cerr << "[ERROR GA] Parent " << i << " is nullptr after selection\n";
                continue;
            }
            for(auto gene : *ind){
                if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                    std::cerr << "[ERROR GA] Parent " << i << " has out-of-range gene after selection: " << gene << "\n";
                }
            }
        }
        
        {
            PROFILE_SCOPE("crossover"); 
            offspring = m_crossover->crossover(parents, instance, m_representation);
        }

        // **Debug Log**: Po crossoverze
        std::cout << "[DEBUG GA] After crossover, offspring size: " << offspring.size() << "\n";
        for(size_t i = 0; i < offspring.size(); i++) {
            auto &ind = offspring[i];
            if(!ind){
                std::cerr << "[ERROR GA] Offspring " << i << " is nullptr after crossover\n";
                continue;
            }
            for(auto gene : *ind){
                if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                    std::cerr << "[ERROR GA] Offspring " << i << " has out-of-range gene after crossover: " << gene << "\n";
                }
            }
        }
        
        {
            PROFILE_SCOPE("mutation");
            m_mutation->mutate(offspring, instance, m_representation);
        }

        // **Debug Log**: Po mutacji
        std::cout << "[DEBUG GA] After mutation, offspring size: " << offspring.size() << "\n";
        for(size_t i = 0; i < offspring.size(); i++) {
            auto &ind = offspring[i];
            if(!ind){
                std::cerr << "[ERROR GA] Offspring " << i << " is nullptr after mutation\n";
                continue;
            }
            for(auto gene : *ind){
                if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                    std::cerr << "[ERROR GA] Offspring " << i << " has out-of-range gene after mutation: " << gene << "\n";
                }
            }
        }

        {
            PROFILE_SCOPE("cache_update");
            m_cache->updatePopulation(offspring, instance, m_fitness, m_representation);
        }
        
        {
            PROFILE_SCOPE("replacement");
            population = m_replacement->replace(population, offspring, instance, m_fitness, m_representation);
        }

        // **Debug Log**: Po zastąpieniu populacji
        std::cout << "[DEBUG GA] After replacement, population size: " << population.size() << "\n";
        for(size_t i = 0; i < population.size(); i++) {
            auto &ind = population[i];
            if(!ind){
                std::cerr << "[ERROR GA] Population " << i << " is nullptr after replacement\n";
                continue;
            }
            for(auto gene : *ind){
                if(gene < 0 || gene >= (int)instance.getSpectrum().size()){
                    std::cerr << "[ERROR GA] Population " << i << " has out-of-range gene after replacement: " << gene << "\n";
                }
            }
        }

        // Szukamy najlepszego osobnika w tym pokoleniu:
        double currentBestFitness = -std::numeric_limits<double>::infinity();
        std::shared_ptr<std::vector<int>> bestInd = nullptr;
        for (auto& ind : population) {
            double fitVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
            if (fitVal > currentBestFitness) {
                currentBestFitness = fitVal;
                bestInd = ind;
            }
        }

        // Jeśli nasz crossover jest adaptacyjny, przekazujemy info o feedbacku
        {
            auto adaptiveCrossover = std::dynamic_pointer_cast<AdaptiveCrossover>(m_crossover);
            if (adaptiveCrossover) {
                adaptiveCrossover->updateFeedback(currentBestFitness);
            }
        }

        // Aktualizacja globalBest co kilka pokoleń
        if (generation % 5 == 0) {
            PROFILE_SCOPE("best");
            updateGlobalBest(population, instance);
        }

        // Możesz też zliczać statystyki w logGenerationStats co np. 10 pokoleń
        // (już w oryginale), ale tutaj zrobimy callback w każdym pokoleniu.

        // Obliczamy coverage i edgeScore dla najlepszego osobnika:
        double coverageVal = 0.0;
        double edgeScoreVal = 0.0;

        // Sprawdzamy, czy nasz IFitness to OptimizedGraphBasedFitness
        auto optFitness = std::dynamic_pointer_cast<OptimizedGraphBasedFitness>(m_fitness);
        if (optFitness && bestInd) {
            const auto& spectrum = instance.getSpectrum();
            int k = instance.getK();
            // Budujemy (lub z cache) graf
            auto graph = optFitness->buildSpectrumGraph(spectrum, k);

            // Tworzymy adjacencyMatrix
            int n = (int)graph.size();
            std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix(
                n, std::vector<PreprocessedEdge>(n)
            );

            for (int i = 0; i < n; i++) {
                for (const auto &edge : graph[i]) {
                    adjacencyMatrix[i][edge.to] = PreprocessedEdge(edge.to, edge.weight, true);
                }
            }

            optFitness->initBuffers(spectrum.size());
            auto analysis = optFitness->analyzePath(*bestInd, adjacencyMatrix);

            // coverage = unikatowe węzły
            coverageVal = analysis.uniqueNodesUsed;
            // edgeScore (w oryginalnym Evaluate to alpha*(edgesWeight1 - 2*edgesWeight2or3) + beta*...) 
            // Tutaj jedynie np. zsumujemy
            // Możesz zmienić w zależności od tego, co chcesz wyświetlać:
            edgeScoreVal = analysis.edgesWeight1 + analysis.edgesWeight2or3;
        }

        // Wywołujemy callback:
        if (progressCallback) {
            progressCallback(generation,
                             m_config.getMaxGenerations(),
                             currentBestFitness,
                             coverageVal,
                             edgeScoreVal,
                             m_theoreticalMaxFitness);
        }

        // **Debug Log**: Najlepszy fitness w bieżącym pokoleniu
        std::cout << "[DEBUG GA] Generation " << generation << ", best fitness: " << currentBestFitness << "\n";

        generation++;
    }

    // Po pętli
    updateGlobalBest(population, instance);
    if (m_globalBestInd) {
        m_bestDNA = m_representation->decodeToDNA(m_globalBestInd, instance);
        m_config.setGlobalBestFitness(m_globalBestFit); // Aktualizacja GAConfig

        // **Debug Log**: Po zakończeniu GA
        std::cout << "[DEBUG GA] Final best fitness: " << m_globalBestFit << ", length of best DNA: " << m_bestDNA.size() << "\n";

        // Jednorazowy callback na koniec (opcjonalnie):
        if (progressCallback) {
            progressCallback(generation,
                             m_config.getMaxGenerations(),
                             m_globalBestFit,
                             0.0, // coverage
                             0.0, // edgeScore
                             m_theoreticalMaxFitness);
        }

        std::cout << "[GA] End after generation " << generation 
                  << ", best-ever fitness = " << m_globalBestFit
                  << ", length of best DNA = " << m_bestDNA.size() << "\n"
                  << "Cache size: " << m_cache->size() << "\n";
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    
    LOG_INFO("Genetic Algorithm completed after " + std::to_string(generation) + 
             " generations (Duration: " + std::to_string(duration) + "s)");
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance &instance)
{
    // Obliczamy teoretyczne maksymalne fitness
    // Dla prostego fitness to po prostu długość DNA
    // Dla bardziej złożonych implementacji można nadpisać tę metodę
    m_theoreticalMaxFitness = instance.getSize();
}

double runGeneticAlgorithmWrapper(const DNAInstance& instance)
{
    GAConfig config;
    if (!config.loadFromFile("config.cfg")) {
        std::cerr << "Failed to load GA configuration\n";
        return 0.0;
    }
    
    // Update instance-specific parameters
    updateConfigWithInstanceParams(instance, config);
    
    // Create and run GA with this config
    auto cache = std::make_shared<CachedPopulation>();
    config.setCache(cache);
    
    GeneticAlgorithm ga(
        config.getRepresentation(),
        config.getSelection(),
        config.getCrossover(config.crossoverType),
        config.getMutation(),
        config.getReplacement(),
        config.getFitness(),
        config.getStopping(),
        cache,
        config
    );
    
    ga.run(instance);
    return ga.getBestFitness();
}

void GeneticAlgorithm::evolve(const DNAInstance& instance) {
    LOG_INFO("Starting evolution");
    LOG_DEBUG("Population size: " + std::to_string(population.size()));
    
    int generation = 0;
    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    std::vector<std::shared_ptr<std::vector<int>>> parents;
    
    while (!m_stopping->stop(population, generation, instance, m_fitness, m_representation)) {
        DEBUG_LOG("Generation " + std::to_string(generation));
        
        {
            PROFILE_SCOPE("selection");
            parents = m_selection->select(population, instance, m_fitness, m_representation);
        }

        {
            PROFILE_SCOPE("crossover");
            offspring = m_crossover->crossover(parents, instance, m_representation);
        }

        {
            PROFILE_SCOPE("mutation");
            m_mutation->mutate(offspring, instance, m_representation);
        }

        {
            PROFILE_SCOPE("cache_update");
            m_cache->updatePopulation(offspring, instance, m_fitness, m_representation);
        }
        
        {
            PROFILE_SCOPE("replacement");
            population = m_replacement->replace(population, offspring, instance, m_fitness, m_representation);
        }

        // Szukamy najlepszego osobnika w tym pokoleniu:
        double currentBestFitness = -std::numeric_limits<double>::infinity();
        std::shared_ptr<std::vector<int>> bestInd = nullptr;
        for (auto& ind : population) {
            double fitVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
            if (fitVal > currentBestFitness) {
                currentBestFitness = fitVal;
                bestInd = ind;
            }
        }

        // Jeśli nasz crossover jest adaptacyjny, przekazujemy info o feedbacku
        {
            auto adaptiveCrossover = std::dynamic_pointer_cast<AdaptiveCrossover>(m_crossover);
            if (adaptiveCrossover) {
                adaptiveCrossover->updateFeedback(currentBestFitness);
            }
        }

        // Aktualizacja globalBest co kilka pokoleń
        if (generation % 5 == 0) {
            PROFILE_SCOPE("best");
            updateGlobalBest(population, instance);
        }

        // Obliczamy coverage i edgeScore dla najlepszego osobnika:
        double coverageVal = 0.0;
        double edgeScoreVal = 0.0;

        // Sprawdzamy, czy nasz IFitness to OptimizedGraphBasedFitness
        auto optFitness = std::dynamic_pointer_cast<OptimizedGraphBasedFitness>(m_fitness);
        if (optFitness && bestInd) {
            const auto& spectrum = instance.getSpectrum();
            int k = instance.getK();
            // Budujemy (lub z cache) graf
            auto graph = optFitness->buildSpectrumGraph(spectrum, k);

            // Tworzymy adjacencyMatrix
            int n = (int)graph.size();
            std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix(
                n, std::vector<PreprocessedEdge>(n)
            );

            for (int i = 0; i < n; i++) {
                for (const auto &edge : graph[i]) {
                    adjacencyMatrix[i][edge.to] = PreprocessedEdge(edge.to, edge.weight, true);
                }
            }

            optFitness->initBuffers(spectrum.size());
            auto analysis = optFitness->analyzePath(*bestInd, adjacencyMatrix);

            // coverage = unikatowe węzły
            coverageVal = analysis.uniqueNodesUsed;
            // edgeScore (w oryginalnym Evaluate to alpha*(edgesWeight1 - 2*edgesWeight2or3) + beta*...) 
            // Tutaj jedynie np. zsumujemy
            // Możesz zmienić w zależności od tego, co chcesz wyświetlać:
            edgeScoreVal = analysis.edgesWeight1 + analysis.edgesWeight2or3;
        }

        // Wywołujemy callback:
        if (progressCallback) {
            progressCallback(generation,
                             m_config.getMaxGenerations(),
                             currentBestFitness,
                             coverageVal,
                             edgeScoreVal,
                             m_theoreticalMaxFitness);
        }

        // **Debug Log**: Najlepszy fitness w bieżącym pokoleniu
        std::cout << "[DEBUG GA] Generation " << generation << ", best fitness: " << currentBestFitness << "\n";

        generation++;
    }
}

double GeneticAlgorithm::getBestFitness() const {
    return m_globalBestFit;
}

std::shared_ptr<std::vector<int>> GeneticAlgorithm::getBestIndividual() const {
    return m_globalBestInd;
}