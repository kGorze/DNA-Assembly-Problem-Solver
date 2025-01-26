#include <gtest/gtest.h>
#include "../include/metaheuristics/representation_impl.h"
#include "../include/dna/dna_instance.h"
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

TEST_F(RepresentationTest, DirectDNARepresentationInitialization_Test) {
    DirectDNARepresentation representation;
    auto instance = std::make_shared<DNAInstance>();
    instance->setK(3);
    instance->setDNA("ACGTACGT");

    auto population = representation.initializePopulation(10, *instance);

    EXPECT_EQ(population.size(), 10);
    for (const auto& individual : population) {
        EXPECT_TRUE(representation.isValid(individual, *instance));
    }
}

TEST_F(RepresentationTest, GraphPathRepresentationInitialization_Test) {
    GraphPathRepresentation representation;
    auto instance = std::make_shared<DNAInstance>();
    instance->setK(3);
    instance->setDNA("ACGTACGT");

    auto population = representation.initializePopulation(10, *instance);

    EXPECT_EQ(population.size(), 10);
    for (const auto& individual : population) {
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

TEST_F(RepresentationTest, EmptyInstance_Test) {
    DirectDNARepresentation directRepresentation;
    GraphPathRepresentation graphRepresentation;
    DNAInstance emptyInstance;

    EXPECT_THROW(directRepresentation.initializePopulation(10, emptyInstance), std::invalid_argument);
    EXPECT_THROW(graphRepresentation.initializePopulation(10, emptyInstance), std::invalid_argument);
}

TEST_F(RepresentationTest, InvalidPopulationSize_Test) {
    DirectDNARepresentation directRepresentation;
    GraphPathRepresentation graphRepresentation;
    auto instance = std::make_shared<DNAInstance>();
    instance->setK(3);
    instance->setDNA("ACGTACGT");

    EXPECT_THROW(directRepresentation.initializePopulation(0, *instance), std::invalid_argument);
    EXPECT_THROW(directRepresentation.initializePopulation(-1, *instance), std::invalid_argument);
    EXPECT_THROW(graphRepresentation.initializePopulation(0, *instance), std::invalid_argument);
    EXPECT_THROW(graphRepresentation.initializePopulation(-1, *instance), std::invalid_argument);
} 