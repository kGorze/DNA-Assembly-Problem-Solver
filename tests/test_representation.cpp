#include <gtest/gtest.h>
#include "metaheuristics/representation_impl.h"
#include "dna/dna_instance.h"
#include <memory>

class RepresentationTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create a test instance
        instance = std::make_unique<DNAInstance>(100, 10, 2, 5, 5, true, 0.8, 0);
        instance->setDNA("ATCGATCGATCGATCGATCG");
        std::vector<std::string> spectrum = {"ATCG", "TCGA", "CGAT", "GATC"};
        instance->setSpectrum(spectrum);
    }
    
    std::unique_ptr<DNAInstance> instance;
};

TEST_F(RepresentationTest, DirectDNARepresentationInitialization) {
    DirectDNARepresentation representation;
    
    // Initialize population
    auto population = representation.initializePopulation(*instance, 10);
    
    // Verify population size
    EXPECT_EQ(population.size(), 10);
    
    // Verify each individual
    for (const auto& individual : population) {
        EXPECT_FALSE(individual->getGenes().empty());
        EXPECT_TRUE(representation.isValid(individual, *instance));
    }
}

TEST_F(RepresentationTest, GraphPathRepresentationInitialization) {
    GraphPathRepresentation representation;
    
    // Initialize population
    auto population = representation.initializePopulation(*instance, 10);
    
    // Verify population size
    EXPECT_EQ(population.size(), 10);
    
    // Verify each individual
    for (const auto& individual : population) {
        EXPECT_FALSE(individual->getGenes().empty());
        EXPECT_TRUE(representation.isValid(individual, *instance));
    }
}

TEST_F(RepresentationTest, DirectDNARepresentationValidity) {
    DirectDNARepresentation representation;
    
    // Create valid individual
    auto individual = std::make_shared<Individual>();
    std::vector<int> validGenes = {0, 1, 2, 3};  // Assuming these indices are valid
    individual->setGenes(validGenes);
    
    // Test validity
    EXPECT_TRUE(representation.isValid(individual, *instance));
    
    // Create invalid individual
    auto invalidIndividual = std::make_shared<Individual>();
    std::vector<int> invalidGenes = {-1, 100, 200};  // Invalid indices
    invalidIndividual->setGenes(invalidGenes);
    
    // Test invalidity
    EXPECT_FALSE(representation.isValid(invalidIndividual, *instance));
}

TEST_F(RepresentationTest, GraphPathRepresentationValidity) {
    GraphPathRepresentation representation;
    
    // Create valid individual
    auto individual = std::make_shared<Individual>();
    std::vector<int> validGenes = {0, 1, 2, 3};  // Assuming these indices form a valid path
    individual->setGenes(validGenes);
    
    // Test validity
    EXPECT_TRUE(representation.isValid(individual, *instance));
    
    // Create invalid individual
    auto invalidIndividual = std::make_shared<Individual>();
    std::vector<int> invalidGenes = {0, 2, 1, 3};  // Assuming these indices don't form a valid path
    invalidIndividual->setGenes(invalidGenes);
    
    // Test invalidity
    EXPECT_FALSE(representation.isValid(invalidIndividual, *instance));
}

TEST_F(RepresentationTest, DirectDNARepresentationToDNA) {
    DirectDNARepresentation representation;
    
    // Create individual
    auto individual = std::make_shared<Individual>();
    std::vector<int> genes = {0, 1, 2, 3};  // Assuming these indices are valid
    individual->setGenes(genes);
    
    // Convert to DNA
    auto dna = representation.toDNA(individual, *instance);
    
    // Verify DNA is not empty
    EXPECT_FALSE(dna.empty());
    
    // Verify DNA contains only valid nucleotides
    for (char c : dna) {
        EXPECT_TRUE(c == 'A' || c == 'T' || c == 'G' || c == 'C');
    }
}

TEST_F(RepresentationTest, GraphPathRepresentationToDNA) {
    GraphPathRepresentation representation;
    
    // Create individual
    auto individual = std::make_shared<Individual>();
    std::vector<int> genes = {0, 1, 2, 3};  // Assuming these indices form a valid path
    individual->setGenes(genes);
    
    // Convert to DNA
    auto dna = representation.toDNA(individual, *instance);
    
    // Verify DNA is not empty
    EXPECT_FALSE(dna.empty());
    
    // Verify DNA contains only valid nucleotides
    for (char c : dna) {
        EXPECT_TRUE(c == 'A' || c == 'T' || c == 'G' || c == 'C');
    }
}

TEST_F(RepresentationTest, EmptyInstance) {
    DirectDNARepresentation directRepresentation;
    GraphPathRepresentation graphRepresentation;
    DNAInstance emptyInstance;
    
    // Test initialization with empty instance
    EXPECT_THROW(directRepresentation.initializePopulation(emptyInstance, 10), std::invalid_argument);
    EXPECT_THROW(graphRepresentation.initializePopulation(emptyInstance, 10), std::invalid_argument);
    
    // Test validity with empty instance
    auto individual = std::make_shared<Individual>();
    EXPECT_FALSE(directRepresentation.isValid(individual, emptyInstance));
    EXPECT_FALSE(graphRepresentation.isValid(individual, emptyInstance));
    
    // Test DNA conversion with empty instance
    EXPECT_TRUE(directRepresentation.toDNA(individual, emptyInstance).empty());
    EXPECT_TRUE(graphRepresentation.toDNA(individual, emptyInstance).empty());
}

TEST_F(RepresentationTest, InvalidPopulationSize) {
    DirectDNARepresentation directRepresentation;
    GraphPathRepresentation graphRepresentation;
    
    // Test initialization with invalid population size
    EXPECT_THROW(directRepresentation.initializePopulation(*instance, 0), std::invalid_argument);
    EXPECT_THROW(directRepresentation.initializePopulation(*instance, -1), std::invalid_argument);
    EXPECT_THROW(graphRepresentation.initializePopulation(*instance, 0), std::invalid_argument);
    EXPECT_THROW(graphRepresentation.initializePopulation(*instance, -1), std::invalid_argument);
} 