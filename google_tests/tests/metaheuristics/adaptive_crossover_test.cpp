#include <gtest/gtest.h>
#include "../base_test.h"
#include "metaheuristics/adaptive_crossover.h"
#include "metaheuristics/crossover_impl.h"
#include "configuration/genetic_algorithm_configuration.h"
#include "dna/dna_instance.h"
#include "metaheuristics/representation.h"
#include "metaheuristics/individual.h"
#include "utils/logging.h"

class AdaptiveCrossoverTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("adaptive_crossover_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();
        // Konfiguracja dla testów
        config = GAConfig();
        // Ustawienie podstawowych parametrów
        config.setPopulationSize(100);
        config.setCrossoverProbability(0.8);
        
        // Tworzenie instancji testowej DNA
        instance = std::make_shared<DNAInstance>();
        instance->setDNA("ACGT");
        instance->setK(2);
        std::vector<std::string> spectrum = {"AC", "CG", "GT"};
        instance->setSpectrum(spectrum);
        
        config.setInstance(instance);
        
        // Inicjalizacja testowanego obiektu
        crossover = std::make_unique<AdaptiveCrossover>(config);
    }

    void TearDown() override {
        BaseTest::TearDown();
        Logger::cleanup();
    }

    GAConfig config;
    std::shared_ptr<DNAInstance> instance;
    std::unique_ptr<AdaptiveCrossover> crossover;
};

// Test inicjalizacji
TEST_F(AdaptiveCrossoverTest, Initialization) {
    ASSERT_NE(crossover, nullptr);
    auto metrics = crossover->getMetrics();
    EXPECT_EQ(metrics.convergenceTime, -1);
    EXPECT_GT(metrics.operatorUsageRates.size(), 0);
}

// Test krzyżowania z poprawnymi rodzicami
TEST_F(AdaptiveCrossoverTest, ValidCrossover) {
    // Przygotowanie rodziców
    std::vector<std::shared_ptr<Individual>> parents;
    std::vector<int> genes1 = {0, 1, 2};
    std::vector<int> genes2 = {2, 1, 0};
    
    auto parent1 = std::make_shared<Individual>(genes1);
    auto parent2 = std::make_shared<Individual>(genes2);
    parents.push_back(parent1);
    parents.push_back(parent2);
    
    // Wykonanie krzyżowania
    auto representation = std::make_shared<PermutationRepresentation>();
    auto offspring = crossover->crossover(parents, *instance, representation);
    
    // Sprawdzenie wyników
    ASSERT_FALSE(offspring.empty());
    for (const auto& child : offspring) {
        ASSERT_NE(child, nullptr);
        EXPECT_EQ(child->getGenes().size(), genes1.size());
    }
}

// Test aktualizacji metryk
TEST_F(AdaptiveCrossoverTest, MetricsUpdate) {
    // Symulacja kilku generacji
    for (int i = 0; i < 5; i++) {
        crossover->updateFeedback(0.5 + i * 0.1); // Rosnące wartości fitness
    }
    
    auto metrics = crossover->getMetrics();
    EXPECT_GT(metrics.bestFitness, 0.0);
    EXPECT_FALSE(metrics.operatorUsageRates.empty());
}
