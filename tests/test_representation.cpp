#include "../include/metaheuristics/representation.h"
#include <gtest/gtest.h>

class RepresentationTest : public ::testing::Test {
protected:
    void SetUp() override {
        instance.setN(10);
        instance.setK(3);
        instance.setDeltaK(1);
        instance.setLNeg(1);
        instance.setLPoz(1);
        instance.setRepAllowed(false);
        instance.setProbablePositive(true);
        instance.setDNA("ACGTACGTAC");
        instance.generateSpectrum();
    }

    DNAInstance instance;
};

TEST_F(RepresentationTest, PermutationRepresentationInitialization_Test) {
    PermutationRepresentation representation;
    auto population = representation.initializePopulation(10, instance);
    ASSERT_EQ(population.size(), 10);
    for (const auto& individual : population) {
        EXPECT_TRUE(representation.isValid(individual, instance));
    }
}

TEST_F(RepresentationTest, PermutationRepresentationValidity) {
    PermutationRepresentation representation;
    auto individual = std::make_shared<Individual>();
    std::vector<int> genes = {0, 1, 2, 3};
    individual->setGenes(genes);
    EXPECT_TRUE(representation.isValid(individual, instance));

    // Test invalid cases
    genes = {0, 1, 1, 3};  // Duplicate value
    individual->setGenes(genes);
    EXPECT_FALSE(representation.isValid(individual, instance));

    genes = {0, 1, 5, 3};  // Out of range value
    individual->setGenes(genes);
    EXPECT_FALSE(representation.isValid(individual, instance));
}

TEST_F(RepresentationTest, PermutationRepresentationToDNA) {
    PermutationRepresentation representation;
    auto individual = std::make_shared<Individual>();
    std::vector<int> genes = {0, 1, 2, 3};
    individual->setGenes(genes);
    
    std::string dna = representation.toDNA(individual, instance);
    EXPECT_FALSE(dna.empty());
    EXPECT_TRUE(dna.find_first_not_of("ACGT") == std::string::npos);
}

TEST_F(RepresentationTest, PopulationInitialization) {
    PermutationRepresentation representation;
    auto population = representation.initializePopulation(10, instance);
    ASSERT_EQ(population.size(), 10);
}

TEST_F(RepresentationTest, IndividualToString) {
    PermutationRepresentation representation;
    auto individual = std::make_shared<Individual>();
    std::vector<int> genes = {0, 1, 2, 3};
    individual->setGenes(genes);
    
    std::string str = representation.toString(individual, instance);
    EXPECT_FALSE(str.empty());
} 