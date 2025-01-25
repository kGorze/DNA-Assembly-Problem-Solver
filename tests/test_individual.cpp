#include <gtest/gtest.h>
#include "metaheuristics/individual.h"
#include <vector>
#include <memory>

TEST(IndividualTest, DefaultConstructor) {
    Individual individual;
    
    // Verify default values
    EXPECT_TRUE(individual.getGenes().empty());
    EXPECT_FALSE(individual.isValid());
    EXPECT_EQ(individual.getFitness(), 0.0);
}

TEST(IndividualTest, SetGenes) {
    Individual individual;
    std::vector<int> genes = {1, 2, 3, 4, 5};
    
    // Set genes
    individual.setGenes(genes);
    
    // Verify genes
    EXPECT_EQ(individual.getGenes(), genes);
}

TEST(IndividualTest, SetFitness) {
    Individual individual;
    double fitness = 0.75;
    
    // Set fitness
    individual.setFitness(fitness);
    
    // Verify fitness
    EXPECT_EQ(individual.getFitness(), fitness);
}

TEST(IndividualTest, Validity) {
    Individual individual;
    
    // Initially invalid
    EXPECT_FALSE(individual.isValid());
    
    // Set valid
    individual.setValid(true);
    EXPECT_TRUE(individual.isValid());
    
    // Set invalid
    individual.setValid(false);
    EXPECT_FALSE(individual.isValid());
}

TEST(IndividualTest, MoveConstructor) {
    Individual original;
    std::vector<int> genes = {1, 2, 3, 4, 5};
    original.setGenes(genes);
    original.setFitness(0.75);
    original.setValid(true);
    
    // Move construct
    Individual moved(std::move(original));
    
    // Verify moved values
    EXPECT_EQ(moved.getGenes(), genes);
    EXPECT_EQ(moved.getFitness(), 0.75);
    EXPECT_TRUE(moved.isValid());
    
    // Original should be in valid but empty state
    EXPECT_TRUE(original.getGenes().empty());
    EXPECT_EQ(original.getFitness(), 0.0);
    EXPECT_FALSE(original.isValid());
}

TEST(IndividualTest, MoveAssignment) {
    Individual original;
    std::vector<int> genes = {1, 2, 3, 4, 5};
    original.setGenes(genes);
    original.setFitness(0.75);
    original.setValid(true);
    
    Individual assigned;
    assigned = std::move(original);
    
    // Verify assigned values
    EXPECT_EQ(assigned.getGenes(), genes);
    EXPECT_EQ(assigned.getFitness(), 0.75);
    EXPECT_TRUE(assigned.isValid());
    
    // Original should be in valid but empty state
    EXPECT_TRUE(original.getGenes().empty());
    EXPECT_EQ(original.getFitness(), 0.0);
    EXPECT_FALSE(original.isValid());
}

TEST(IndividualTest, InvalidFitness) {
    Individual individual;
    
    // Test setting invalid fitness values
    EXPECT_THROW(individual.setFitness(-1.0), std::invalid_argument);
    EXPECT_THROW(individual.setFitness(1.1), std::invalid_argument);
    EXPECT_THROW(individual.setFitness(std::numeric_limits<double>::infinity()), std::invalid_argument);
    EXPECT_THROW(individual.setFitness(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

TEST(IndividualTest, EmptyGenes) {
    Individual individual;
    std::vector<int> emptyGenes;
    
    // Set empty genes
    individual.setGenes(emptyGenes);
    
    // Verify empty genes
    EXPECT_TRUE(individual.getGenes().empty());
    EXPECT_FALSE(individual.isValid());
}

TEST(IndividualTest, CopyDisabled) {
    Individual individual;
    std::vector<int> genes = {1, 2, 3, 4, 5};
    individual.setGenes(genes);
    
    // These lines should not compile
    // Individual copy(individual);  // Copy constructor
    // Individual assigned = individual;  // Copy assignment
    
    // Verify that we can still move
    Individual moved(std::move(individual));
    EXPECT_EQ(moved.getGenes(), genes);
} 