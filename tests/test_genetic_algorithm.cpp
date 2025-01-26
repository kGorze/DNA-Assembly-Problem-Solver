#include <gtest/gtest.h>
#include "../include/metaheuristics/genetic_algorithm.h"
#include "../include/dna/dna_instance.h"
#include "../include/metaheuristics/representation_impl.h"
#include <memory>

class GeneticAlgorithmTest : public ::testing::Test {
protected:
    void SetUp() override {
        instance.setK(4);
        instance.setDNA("ACGTACGT");
        instance.setSpectrum({"ACGT", "CGTA", "GTAC", "TACG"});
        
        config.setPopulationSize(50);
        config.setMaxGenerations(100);
        config.setMutationRate(0.1);
        config.setCrossoverProbability(0.8);
        config.setTournamentSize(5);
        config.setTargetFitness(1.0);
        
        representation = std::make_unique<DirectDNARepresentation>();
    }
    
    DNAInstance instance;
    GAConfig config;
    std::unique_ptr<IRepresentation> representation;
};

TEST_F(GeneticAlgorithmTest, InitializePopulation) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    std::string result = algorithm.run(instance);
    EXPECT_FALSE(result.empty());
}

TEST_F(GeneticAlgorithmTest, RunMultipleGenerations) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    std::string result = algorithm.run(instance);
    EXPECT_FALSE(result.empty());
    EXPECT_GT(algorithm.getBestFitness(), 0.0);
}

TEST_F(GeneticAlgorithmTest, TargetFitnessTermination) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    std::string result = algorithm.run(instance);
    EXPECT_FALSE(result.empty());
    EXPECT_LE(algorithm.getBestFitness(), config.getTargetFitness());
}

TEST_F(GeneticAlgorithmTest, InvalidConfig) {
    GAConfig invalidConfig;
    invalidConfig.setPopulationSize(0);  // Invalid population size
    invalidConfig.setMaxGenerations(0);  // Invalid generations
    invalidConfig.setMutationRate(-0.1); // Invalid mutation rate
    
    EXPECT_THROW({
        GeneticAlgorithm algorithm(std::move(representation), invalidConfig);
    }, std::invalid_argument);
}

TEST_F(GeneticAlgorithmTest, EmptyInstance) {
    DNAInstance emptyInstance;
    GeneticAlgorithm algorithm(std::move(representation), config);
    std::string result = algorithm.run(emptyInstance);
    EXPECT_TRUE(result.empty());
}

TEST_F(GeneticAlgorithmTest, GraphPathRepresentation) {
    auto graphRepresentation = std::make_unique<GraphPathRepresentation>();
    GeneticAlgorithm algorithm(std::move(graphRepresentation), config);
    std::string result = algorithm.run(instance);
    EXPECT_FALSE(result.empty());
}

TEST_F(GeneticAlgorithmTest, ConfigValidation) {
    GAConfig testConfig = config;
    
    // Test invalid population size
    testConfig.setPopulationSize(0);
    EXPECT_THROW({
        GeneticAlgorithm algorithm(std::move(representation), testConfig);
    }, std::invalid_argument);
    
    // Reset config
    testConfig = config;
    
    // Test invalid mutation rate
    testConfig.setMutationRate(-0.1);
    EXPECT_THROW({
        GeneticAlgorithm algorithm(std::move(representation), testConfig);
    }, std::invalid_argument);
    
    // Reset config
    testConfig = config;
    
    // Test invalid crossover probability
    testConfig.setCrossoverProbability(1.5);
    EXPECT_THROW({
        GeneticAlgorithm algorithm(std::move(representation), testConfig);
    }, std::invalid_argument);
}

TEST_F(GeneticAlgorithmTest, ResultValidation) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    std::string result = algorithm.run(instance);
    
    // Check if result is a valid DNA sequence
    EXPECT_FALSE(result.empty());
    for (char c : result) {
        EXPECT_TRUE(c == 'A' || c == 'C' || c == 'G' || c == 'T');
    }
    
    // Check if length is reasonable
    EXPECT_GE(result.length(), instance.getK());
    
    // Check if fitness is within bounds
    double fitness = algorithm.getBestFitness();
    EXPECT_GE(fitness, 0.0);
    EXPECT_LE(fitness, 1.0);
} 