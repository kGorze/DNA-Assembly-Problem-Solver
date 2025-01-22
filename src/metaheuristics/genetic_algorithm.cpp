#include "metaheuristics/genetic_algorithm.h"
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include <mutex>

std::mutex GeneticAlgorithm::outputMutex;


GeneticAlgorithm::GeneticAlgorithm(
    std::shared_ptr<IRepresentation> representation,
    std::shared_ptr<ISelection> selection,
    std::shared_ptr<ICrossover> crossover,
    std::shared_ptr<IMutation> mutation,
    std::shared_ptr<IReplacement> replacement,
    std::shared_ptr<IFitness> fitness,
    std::shared_ptr<IStopping> stopping,
    std::shared_ptr<IPopulationCache> cache)
    : m_representation(representation)
    , m_selection(selection)
    , m_crossover(crossover)
    , m_mutation(mutation)
    , m_replacement(replacement)
    , m_fitness(fitness)
    , m_stopping(stopping)
    , m_cache(cache)
    , m_globalBestInd(nullptr)
    , m_globalBestFit(-std::numeric_limits<double>::infinity())
{ }

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
    double progress = (generation * 100.0) / GAConfig::getInstance().getMaxGenerations();
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
}

void GeneticAlgorithm::updateGlobalBest(
    const std::vector<std::shared_ptr<std::vector<int>>> &pop,
    const DNAInstance &instance)
{
    for (const auto &ind : pop) {
        if (!ind) continue;
        double fitVal = m_cache->getOrCalculateFitness(ind, instance, m_fitness, m_representation);
        if (fitVal > m_globalBestFit) {
            m_globalBestFit = fitVal;
            m_globalBestInd = std::make_shared<std::vector<int>>(*ind);
        }
    }
}

void GeneticAlgorithm::run(const DNAInstance &instance)
{
    PROFILE_FUNCTION();

    population = m_representation->initializePopulation(GAConfig::getInstance().getPopulationSize(), instance);
    m_cache->updatePopulation(population, instance, m_fitness, m_representation);

    int generation = 0;
    int maxGenerations = GAConfig::getInstance().getMaxGenerations();

    // Z góry obliczamy teoretyczne maksymalne fitness (np. w oparciu o założenia w OptimizedGraphBasedFitness)
    calculateTheoreticalMaxFitness(instance);

    std::vector<std::shared_ptr<std::vector<int>>> offspring;
    std::vector<std::shared_ptr<std::vector<int>>> parents;
    
    while (!m_stopping->stop(population, generation, instance, m_fitness, m_representation)) {
        
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

        // Możemy też zliczać statystyki w logGenerationStats co np. 10 pokoleń
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
            std::vector<std::vector<OptimizedGraphBasedFitness::PreprocessedEdge>> adjacencyMatrix(
                n, std::vector<OptimizedGraphBasedFitness::PreprocessedEdge>(n)
            );

            for (int i = 0; i < n; i++) {
                for (auto &edge : graph[i]) {
                    adjacencyMatrix[i][edge.to] = { edge.to, edge.weight, true };
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
                             maxGenerations,
                             currentBestFitness,
                             coverageVal,
                             edgeScoreVal,
                             m_theoreticalMaxFitness);
        }

        generation++;
    }

    // Po zakończeniu pętli stop:
    updateGlobalBest(population, instance);
    if (m_globalBestInd) {
        m_bestDNA = m_representation->decodeToDNA(m_globalBestInd, instance);

        // Jednorazowy callback na koniec (opcjonalnie):
        if (progressCallback) {
            progressCallback(generation,
                             maxGenerations,
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
}

void GeneticAlgorithm::calculateTheoreticalMaxFitness(const DNAInstance &instance) {
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();

    if (spectrum.empty() || k <= 0) {
        m_theoreticalMaxFitness = 0.0;
        return;
    }

    size_t n = spectrum.size();

    // Zakładamy, że "best edge score" to (n-1), a coverage to n
    double bestEdgeScore = (n - 1);
    double bestCoverageScore = n;

    // Współczynniki domyślne w OptimizedGraphBasedFitness:
    double alpha = 0.7; 
    double beta = 0.3;  

    m_theoreticalMaxFitness = (alpha * bestEdgeScore) + (beta * bestCoverageScore);

    std::cout << "[MaxFitness] Theoretical maximum calculated:\n"
              << "- Best edge score: " << bestEdgeScore << " (weight: " << alpha << ")\n"
              << "- Best coverage score: " << bestCoverageScore << " (weight: " << beta << ")\n"
              << "- Total theoretical maximum: " << m_theoreticalMaxFitness << "\n";
}
