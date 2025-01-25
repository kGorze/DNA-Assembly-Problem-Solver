#include <gtest/gtest.h>
#include "metaheuristics/genetic_algorithm.h"
#include "dna/dna_instance.h"
#include "metaheuristics/representation_impl.h"
#include <memory>

class GeneticAlgorithmTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a test instance
        instance = std::make_unique<DNAInstance>(100, 10, 2, 5, 5, true, 0.8, 0);
        instance->setDNA("ATCGATCGATCGATCGATCG");
        std::vector<std::string> spectrum = {"ATCG", "TCGA", "CGAT", "GATC"};
        instance->setSpectrum(spectrum);
        
        // Create genetic config
        config.populationSize = 50;
        config.maxGenerations = 100;
        config.mutationProbability = 0.1;
        config.crossoverProbability = 0.8;
        config.tournamentSize = 3;
        config.targetFitness = 0.95;
        
        // Create representation
        representation = std::make_unique<DirectDNARepresentation>();
    }
    
    std::unique_ptr<DNAInstance> instance;
    GeneticConfig config;
    std::unique_ptr<IRepresentation> representation;
};

TEST_F(GeneticAlgorithmTest, InitializePopulation) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    
    // Run algorithm for 1 generation
    config.maxGenerations = 1;
    auto result = algorithm.run(*instance);
    
    // Verify result is not empty
    EXPECT_FALSE(result.empty());
}

TEST_F(GeneticAlgorithmTest, RunMultipleGenerations) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    
    // Run algorithm for multiple generations
    config.maxGenerations = 10;
    auto result = algorithm.run(*instance);
    
    // Verify result is not empty
    EXPECT_FALSE(result.empty());
}

TEST_F(GeneticAlgorithmTest, TargetFitnessTermination) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    
    // Set a low target fitness to ensure early termination
    config.targetFitness = 0.1;
    config.maxGenerations = 1000;
    
    auto result = algorithm.run(*instance);
    
    // Verify result is not empty
    EXPECT_FALSE(result.empty());
}

TEST_F(GeneticAlgorithmTest, InvalidConfig) {
    // Create invalid config
    GeneticConfig invalidConfig;
    invalidConfig.populationSize = 0;  // Invalid population size
    invalidConfig.maxGenerations = -1;  // Invalid generation count
    invalidConfig.mutationProbability = 1.5;  // Invalid probability
    invalidConfig.crossoverProbability = -0.1;  // Invalid probability
    invalidConfig.tournamentSize = 0;  // Invalid tournament size
    
    EXPECT_THROW(GeneticAlgorithm(std::move(representation), invalidConfig), std::invalid_argument);
}

TEST_F(GeneticAlgorithmTest, EmptyInstance) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    
    // Create empty instance
    DNAInstance emptyInstance;
    
    // Run algorithm with empty instance
    EXPECT_THROW(algorithm.run(emptyInstance), std::invalid_argument);
}

TEST_F(GeneticAlgorithmTest, GraphPathRepresentation) {
    // Create a new instance of GraphPathRepresentation
    auto graphRepresentation = std::make_unique<GraphPathRepresentation>();
    
    GeneticAlgorithm algorithm(std::move(graphRepresentation), config);
    
    // Run algorithm
    auto result = algorithm.run(*instance);
    
    // Verify result is not empty
    EXPECT_FALSE(result.empty());
}

TEST_F(GeneticAlgorithmTest, ConfigValidation) {
    // Test various invalid configurations
    GeneticConfig testConfig = config;
    
    // Test population size
    testConfig.populationSize = 0;
    EXPECT_THROW(GeneticAlgorithm(std::move(std::make_unique<DirectDNARepresentation>()), testConfig), std::invalid_argument);
    testConfig = config;
    
    // Test max generations
    testConfig.maxGenerations = 0;
    EXPECT_THROW(GeneticAlgorithm(std::move(std::make_unique<DirectDNARepresentation>()), testConfig), std::invalid_argument);
    testConfig = config;
    
    // Test mutation probability
    testConfig.mutationProbability = -0.1;
    EXPECT_THROW(GeneticAlgorithm(std::move(std::make_unique<DirectDNARepresentation>()), testConfig), std::invalid_argument);
    testConfig = config;
    
    // Test crossover probability
    testConfig.crossoverProbability = 1.1;
    EXPECT_THROW(GeneticAlgorithm(std::move(std::make_unique<DirectDNARepresentation>()), testConfig), std::invalid_argument);
    testConfig = config;
    
    // Test tournament size
    testConfig.tournamentSize = 0;
    EXPECT_THROW(GeneticAlgorithm(std::move(std::make_unique<DirectDNARepresentation>()), testConfig), std::invalid_argument);
}

TEST_F(GeneticAlgorithmTest, ResultValidation) {
    GeneticAlgorithm algorithm(std::move(representation), config);
    
    // Run algorithm
    auto result = algorithm.run(*instance);
    
    // Verify result format
    EXPECT_TRUE(result.find("Generation") != std::string::npos);
    EXPECT_TRUE(result.find("Best Fitness") != std::string::npos);
    EXPECT_TRUE(result.find("Average Fitness") != std::string::npos);
} 