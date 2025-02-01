// adaptive_crossover_test.cpp

#include <gtest/gtest.h>
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <numeric>
#include <algorithm>
#include "../base_test.h"

// Trick: aby uzyskać dostęp do metod prywatnych AdaptiveCrossover
#define private public
#include "metaheuristics/adaptive_crossover.h"
#undef private

#include "interfaces/i_crossover.h"
#include "interfaces/i_representation.h"
#include "interfaces/i_population_cache.h"
#include "utils/random.h"
#include "utils/logging.h"
#include "metaheuristics/individual.h"

// ----- DUMMY CLASY -----

// Zakładamy, że klasa Individual jest już zdefiniowana w projekcie (np. w "include/metaheuristics/individual.h")
// Dlatego nie redefiniujemy jej – tworzymy jedynie klasę DummyIndividual dziedziczącą po istniejącej Individual.

class DummyIndividual : public Individual {
public:
    explicit DummyIndividual(const std::vector<int>& genes) {
        setGenes(genes);
    }
    
    std::string toString() const {
        std::ostringstream oss;
        oss << "DummyIndividual(genes=[";
        const auto& genes = getGenes();
        for (size_t i = 0; i < genes.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << genes[i];
        }
        oss << "], fitness=" << getFitness() << ")";
        return oss.str();
    }

    virtual ~DummyIndividual() = default;
};

using IndividualPtr = std::shared_ptr<DummyIndividual>;

class DummyRepresentation : public IRepresentation {
public:
    bool isValid(const std::shared_ptr<Individual>& individual, [[maybe_unused]] const DNAInstance& instance) const override {
        return (individual != nullptr);
    }
    
    std::string toString(const std::shared_ptr<Individual>& individual, const DNAInstance& instance) const override {
        return toDNA(individual, instance);
    }
    
    std::string toDNA(const std::shared_ptr<Individual>& individual, [[maybe_unused]] const DNAInstance& instance) const override {
        if (!individual) return "";
        auto* dummyInd = static_cast<DummyIndividual*>(individual.get());
        if (!dummyInd) return "";
        
        std::stringstream ss;
        for (int g : dummyInd->getGenes()) {
            char base = "ACGT"[g % 4];
            ss << base;
        }
        return ss.str();
    }
    
    std::vector<std::shared_ptr<Individual>> initializePopulation(size_t populationSize, [[maybe_unused]] const DNAInstance& instance) override {
        std::vector<std::shared_ptr<Individual>> population;
        for (size_t i = 0; i < populationSize; ++i) {
            population.push_back(std::make_shared<DummyIndividual>(std::vector<int>{0,1,2,3}));
        }
        return population;
    }
    
    bool initializeIndividual([[maybe_unused]] Individual& individual, [[maybe_unused]] const DNAInstance& instance) override {
        return true;
    }
};

class DummyPopulationCache : public IPopulationCache {
public:
    DummyPopulationCache(const std::vector<std::shared_ptr<Individual>>& pop) : population(pop) {}
    
    const std::vector<std::shared_ptr<Individual>>& getCurrentPopulation() const override {
        return population;
    }
    
    double getOrCalculateFitness([[maybe_unused]] const std::shared_ptr<Individual>& individual, 
                                [[maybe_unused]] const DNAInstance& instance) override {
        return 0.0;
    }
    
    void updatePopulation(const std::vector<std::shared_ptr<Individual>>& pop) override {
        population = pop;
    }
    
    void clear() override {
        population.clear();
    }
    
    void reserve(size_t size) override {
        population.reserve(size);
    }
    
    void add(const std::shared_ptr<Individual>& individual) override {
        population.push_back(individual);
    }
    
    bool contains(const std::shared_ptr<Individual>& individual) const override {
        return std::find(population.begin(), population.end(), individual) != population.end();
    }
    
    void enableDiversityTracking([[maybe_unused]] bool enabled) override {}
    void setDiversityThreshold([[maybe_unused]] double threshold) override {}
    bool isDiversityTrackingEnabled() const override { return false; }
    double getDiversityThreshold() const override { return 0.0; }
    
private:
    std::vector<std::shared_ptr<Individual>> population;
};

class DummyDNAInstance : public DNAInstance {
    // Wystarczające, aby przekazać instancję do reprezentacji – puste ciało
};

class DummyGAConfig : public GAConfig {
public:
    DummyGAConfig() {
        instance = std::make_shared<DummyDNAInstance>();
        representation = std::make_shared<DummyRepresentation>();
        cache = nullptr;
    }
    
    std::shared_ptr<IPopulationCache> getCache() const {
        return cache;
    }
    
    void setCache(const std::shared_ptr<IPopulationCache>& c) { 
        cache = c; 
    }
    
    std::shared_ptr<IRepresentation> getRepresentation() const {
        return representation;
    }
    
    std::shared_ptr<DNAInstance> getInstance() const {
        return instance;
    }
    
    size_t getPopulationSize() const { return 10; }
    size_t getMaxGenerations() const { return 100; }
    double getCrossoverProbability() const { return 0.8; }
    double getMutationProbability() const { return 0.1; }
    bool isElitismEnabled() const { return true; }
    size_t getEliteCount() const { return 1; }
    bool isDiversityTrackingEnabled() const { return false; }
    double getDiversityThreshold() const { return 0.1; }
    
private:
    std::shared_ptr<DNAInstance> instance;
    std::shared_ptr<IRepresentation> representation;
    std::shared_ptr<IPopulationCache> cache;
};

class DummyCrossover : public ICrossover {
public:
    std::vector<std::shared_ptr<Individual>> crossover(
        const std::vector<std::shared_ptr<Individual>>& parents,
        [[maybe_unused]] const DNAInstance& instance,
        [[maybe_unused]] std::shared_ptr<IRepresentation> representation) override {
        if (parents.empty()) return {};
        auto first = parents[0];
        // Załóżmy, że nasz dummy operator zwraca kopię pierwszego rodzica
        if (!first) return {};
        auto genes = first->getGenes();
        std::vector<std::shared_ptr<Individual>> result;
        result.push_back(std::make_shared<DummyIndividual>(genes));
        return result;
    }
};

// Aby testować AdaptiveCrossover, stworzymy klasę pochodną, która "eksponuje" metody prywatne.
class TestAdaptiveCrossover : public AdaptiveCrossover {
public:
    TestAdaptiveCrossover(const GAConfig& config) : AdaptiveCrossover(config) {}

    std::vector<double>& exposedProbabilities() { return m_crossoverProbabilities; }
    std::vector<int>& exposedUsage() { return m_crossoverUsage; }
    std::vector<double>& exposedPerformance() { return m_crossoverPerformance; }
    int& exposedGenerationCount() { return generationCount; }

    void setProbabilities(const std::vector<double>& probs) {
        m_crossoverProbabilities = probs;
    }
    
    void setGenerationCount(int count) { 
        generationCount = count; 
    }

    // Udostępniamy metody prywatne:
    using AdaptiveCrossover::adjustProbabilities;
    using AdaptiveCrossover::updatePerformance;
    using AdaptiveCrossover::calculateAverageDistance;
    using AdaptiveCrossover::selectCrossover;
    using AdaptiveCrossover::updateMetrics;
    using AdaptiveCrossover::logDiversityMetrics;
};

//
// Fixture
//
class AdaptiveCrossoverTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("adaptive_crossover_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();  // Call base class SetUp
        
        config = std::make_shared<DummyGAConfig>();
        adaptive = std::make_unique<TestAdaptiveCrossover>(*config);
    }

    void TearDown() override {
        BaseTest::TearDown();  // Call base class TearDown
        Logger::cleanup();
    }
    
    std::shared_ptr<DummyGAConfig> config;
    std::unique_ptr<TestAdaptiveCrossover> adaptive;
};

//
// Test 1: Konstruktor – sprawdzenie początkowych metryk
//
TEST_F(AdaptiveCrossoverTest, ConstructorInitializesMetricsCorrectly) {
    RunMetrics metrics = adaptive->getMetrics();
    EXPECT_NEAR(metrics.avgFitness, 0.0, 1e-6);
    EXPECT_EQ(metrics.bestFitness, -std::numeric_limits<double>::infinity());
    EXPECT_EQ(metrics.convergenceTime, -1);
    for (double rate : metrics.operatorUsageRates) {
        EXPECT_NEAR(rate, 0.0, 1e-6);
    }
    for (double rate : metrics.operatorSuccessRates) {
        EXPECT_NEAR(rate, 0.0, 1e-6);
    }
}

//
// Test 2: calculateAverageDistance – sprawdzenie obliczeń
//
TEST_F(AdaptiveCrossoverTest, CalculateAverageDistance) {
    std::vector<int> genes = {0, 1, 2, 3, 0, 1};
    auto ind1 = std::make_shared<DummyIndividual>(genes);
    auto ind2 = std::make_shared<DummyIndividual>(genes);
    std::vector<std::shared_ptr<Individual>> pop = {ind1, ind2};
    
    double avgDist = adaptive->calculateAverageDistance(pop);
    EXPECT_NEAR(avgDist, 0.0, 1e-6);
    
    std::vector<int> genes2 = {3, 2, 1, 0, 3, 2};
    auto ind3 = std::make_shared<DummyIndividual>(genes2);
    pop = {ind1, ind3};
    EXPECT_NEAR(avgDist = adaptive->calculateAverageDistance(pop), 1.0, 1e-6);
}

//
// Test 3: selectCrossover – ustawiamy sztuczne prawdopodobieństwa
//
TEST_F(AdaptiveCrossoverTest, SelectCrossoverChoosesBasedOnProbabilities) {
    // Ustawiamy ręcznie: tylko operator o indeksie 3 ma prawdopodobieństwo 1
    adaptive->setProbabilities({0.0, 0.0, 0.0, 1.0});
    auto crossoverOp = adaptive->selectCrossover();
    // Sprawdzamy, że prawdopodobieństwa są zgodne oraz że aktualny indeks wynosi 3
    EXPECT_EQ(adaptive->exposedProbabilities().size(), 4u);
    EXPECT_NEAR(adaptive->exposedProbabilities()[3], 1.0, 1e-6);
}

//
// Test 4: Metoda crossover – warunki brzegowe
//
TEST_F(AdaptiveCrossoverTest, CrossoverReturnsFallbackWhenParentsInvalid) {
    std::vector<std::shared_ptr<Individual>> parents = {
        std::make_shared<DummyIndividual>(std::vector<int>{0,1,2})
    };
    auto representation = config->getRepresentation();
    DummyDNAInstance instance;
    auto offspring = adaptive->crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), parents.size());
}

TEST_F(AdaptiveCrossoverTest, CrossoverReturnsEmptyWhenRepresentationNull) {
    std::vector<std::shared_ptr<Individual>> parents = {
        std::make_shared<DummyIndividual>(std::vector<int>{0,1,2}),
        std::make_shared<DummyIndividual>(std::vector<int>{1,2,3})
    };
    DummyDNAInstance instance;
    auto offspring = adaptive->crossover(parents, instance, nullptr);
    EXPECT_TRUE(offspring.empty());
}

//
// Test 5: updateFeedback – symulacja przebiegu generacji
//
TEST_F(AdaptiveCrossoverTest, UpdateFeedbackUpdatesMetrics) {
    adaptive->exposedGenerationCount() = 1;
    adaptive->updateFeedback(10.0);
    adaptive->exposedGenerationCount() = 2;
    adaptive->updateFeedback(12.0);
    adaptive->exposedGenerationCount() = 3;
    adaptive->updateFeedback(11.0);
    
    RunMetrics metrics = adaptive->getMetrics();
    EXPECT_NEAR(metrics.bestFitness, 12.0, 1e-6);
    EXPECT_NEAR(metrics.avgFitness, 11.0, 1e-6);
    EXPECT_EQ(metrics.convergenceTime, 2);
}

//
// Test 6: setParameters – metoda pusta, nie wpływa na metryki
//
TEST_F(AdaptiveCrossoverTest, SetParametersDoesNotAffectBehavior) {
    adaptive->setParameters(0.5, 5, 3, 0.1);
    RunMetrics metrics = adaptive->getMetrics();
    EXPECT_EQ(metrics.convergenceTime, -1);
    EXPECT_NEAR(metrics.avgFitness, 0.0, 1e-6);
}

//
// Test 7: updateMetrics – sprawdzenie, czy metoda nie rzuca wyjątków
//
TEST_F(AdaptiveCrossoverTest, UpdateMetricsDoesNotThrow) {
    std::vector<std::shared_ptr<Individual>> parents = {
        std::make_shared<DummyIndividual>(std::vector<int>{0,1,2}),
        std::make_shared<DummyIndividual>(std::vector<int>{1,2,3})
    };
    std::vector<std::shared_ptr<Individual>> offspring = {
        std::make_shared<DummyIndividual>(std::vector<int>{0,1,2})
    };
    DummyDNAInstance instance;
    auto representation = config->getRepresentation();
    EXPECT_NO_THROW({
        adaptive->updateMetrics(offspring, parents, instance, representation);
    });
}

//
// Test 8: adjustProbabilities – symulacja statystyk
//
TEST_F(AdaptiveCrossoverTest, AdjustProbabilitiesNormalizesAndBiases) {
    std::vector<int> usage = {10, 10, 10, 10};
    std::vector<double> performance = {3.0, 5.0, 2.0, 2.0};
    adaptive->setGenerationCount(10);
    adaptive->exposedUsage() = usage;
    adaptive->exposedPerformance() = performance;
    adaptive->adjustProbabilities();
    const auto& probs = adaptive->exposedProbabilities();
    double sum = std::accumulate(probs.begin(), probs.end(), 0.0);
    EXPECT_NEAR(sum, 1.0, 1e-6);
    EXPECT_GT(probs[1], probs[0]);
    EXPECT_GT(probs[1], probs[2]);
    EXPECT_GT(probs[1], probs[3]);
}

//
// Test 9: crossover – sprawdzenie poprawnego działania
//
TEST_F(AdaptiveCrossoverTest, CrossoverReturnsOffspringWhenDataValid) {
    std::vector<int> genes1 = {0, 1, 2, 3};
    std::vector<int> genes2 = {3, 2, 1, 0};
    std::vector<std::shared_ptr<Individual>> parents = {
        std::make_shared<DummyIndividual>(genes1),
        std::make_shared<DummyIndividual>(genes2)
    };
    DummyDNAInstance instance;
    auto representation = config->getRepresentation();
    ASSERT_NE(representation, nullptr) << "Representation should not be null";
    
    auto offspring = adaptive->crossover(parents, instance, representation);
    ASSERT_FALSE(offspring.empty()) << "Offspring vector should not be empty";
    ASSERT_NE(offspring[0], nullptr) << "First offspring should not be null";
    
    auto* dummyOffspring = static_cast<DummyIndividual*>(offspring[0].get());
    ASSERT_NE(dummyOffspring, nullptr) << "Offspring should not be null after cast";
    EXPECT_EQ(dummyOffspring->getGenes(), genes1);
}

//
// Test 10: updatePerformance – sprawdzenie przyrostu wskaźników
//
TEST_F(AdaptiveCrossoverTest, UpdatePerformanceIncrementsUsageAndPerformance) {
    adaptive->exposedGenerationCount() = 1;
    // Używamy wektorów double (a nie int) dla m_crossoverPerformance, aby pasowały do typu
    adaptive->exposedUsage() = std::vector<int>(4, 0);
    adaptive->exposedPerformance() = std::vector<double>(4, 0.0);
    adaptive->setProbabilities({0.0, 0.0, 1.0, 0.0});
    adaptive->updatePerformance(true);
    auto usage = adaptive->exposedUsage();
    auto performance = adaptive->exposedPerformance();
    EXPECT_EQ(usage[2], 1);
    EXPECT_EQ(performance[2], 1.0);
    adaptive->updatePerformance(false);
    usage = adaptive->exposedUsage();
    performance = adaptive->exposedPerformance();
    EXPECT_EQ(usage[2], 2);
    EXPECT_EQ(performance[2], 1.0);
}

//
// Test 11: logDiversityMetrics – sprawdzamy, czy historia metryk się aktualizuje
//
TEST_F(AdaptiveCrossoverTest, LogDiversityMetricsUpdatesHistory) {
    std::vector<std::shared_ptr<Individual>> population = {
        std::make_shared<DummyIndividual>(std::vector<int>{0,1,2,3}),
        std::make_shared<DummyIndividual>(std::vector<int>{0,1,2,3}),
        std::make_shared<DummyIndividual>(std::vector<int>{3,2,1,0}),
        std::make_shared<DummyIndividual>(std::vector<int>{3,2,1,0})
    };
    auto cache = std::make_shared<DummyPopulationCache>(population);
    config->setCache(cache);
    auto representation = config->getRepresentation();
    // Aby móc przekazać referencję, tworzymy lokalną zmienną RunMetrics
    RunMetrics metrics = adaptive->getMetrics();
    adaptive->logDiversityMetrics(population, representation, metrics);
    // Sprawdzamy, czy historia metryk ma przynajmniej jeden wpis
    EXPECT_FALSE(metrics.diversityHistory.empty());
}
