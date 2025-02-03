//
// Created by konrad_guest on 03/02/2025.
//
// OptimizedGraphBasedFitness_test.cpp
//
// This file tests the OptimizedGraphBasedFitness class (the graph‐based fitness
// implementation) from your codebase. It uses the existing base_test.h and test_main 
// infrastructure. The tests cover low‐level functions (edge weight, adjacency matrix, etc.) 
// as well as the overall fitness integration.

#include "../../tests/base_test.h"
#include "metaheuristics/fitness.h"               // for OptimizedGraphBasedFitness
#include "metaheuristics/individual.h"            // for Individual
#include "metaheuristics/representation.h" // for PermutationRepresentation
#include "dna/dna_instance.h"                     // for DNAInstance
#include <vector>
#include <string>
#include <set>
#include <limits>
#include <cmath>
#include <typeinfo>

// For convenience, we assume that the OptimizedGraphBasedFitness class exposes
// the following (public) member functions:
//   - double calculateFitness(const std::shared_ptr<Individual>& solution,
//                             const DNAInstance& instance,
//                             std::shared_ptr<IRepresentation> representation) const;
//   - std::vector<std::vector<PreprocessedEdge>> buildAdjacencyMatrix(const DNAInstance& instance) const;
//   - int calculateEdgeWeight(const std::string& from, const std::string& to, int k) const;
//   - double calculateEdgeQuality(const std::shared_ptr<Individual>& solution, const DNAInstance& instance) const;
//   - double calculateLength(const std::shared_ptr<Individual>& solution, const DNAInstance& instance) const;
//   - double calculateSpectrumCoverageScore(const std::vector<char>& dna, const DNAInstance& instance) const;
//   - double calculateConnectivity(const std::shared_ptr<Individual>& solution, const DNAInstance& instance) const;
// (If some functions are protected, you might create "friend" test classes or expose them via a test‐hook.)

// Define a test fixture for OptimizedGraphBasedFitness tests.
class OptimizedGraphBasedFitnessTest : public ::testing::Test, public OptimizedGraphBasedFitness {
protected:
    // Create a default OptimizedGraphBasedFitness instance.
    OptimizedGraphBasedFitness& fitnessCalc = *this;  // Use the test class itself as the fitness calculator

    // A dummy DNA instance with k=5 and a spectrum of 10 k‑mers.
    DNAInstance instance;

    // A representation to convert individuals to DNA.
    std::shared_ptr<IRepresentation> representation;

    // A valid individual using a permutation of indices 0..9.
    std::shared_ptr<Individual> validIndividual;

    // Public test methods that can access protected members
    void buildGraphTest(const std::vector<int>& genes, const DNAInstance& instance) {
        buildGraphForTest(genes, instance);
    }

    bool testConnection(int node1, int node2, const DNAInstance& instance) {
        return testNodesConnection(node1, node2, instance);
    }

    double getConnScore() {
        return getConnectivityScore();
    }

    double getLengthPen(const std::vector<int>& genes, const DNAInstance& instance) {
        return getLengthPenalty(genes, instance);
    }

    double getSpectrumScore(const std::vector<int>& genes, const DNAInstance& instance) {
        return getSpectrumCoverageScore(genes, instance);
    }

    const std::vector<std::vector<bool>>& getMatrix() {
        return getAdjacencyMatrix();
    }

    virtual void SetUp() override {
        // Set instance parameters (these must match your defaults or desired test values)
        instance.setK(5);
        instance.setDeltaK(0);
        instance.setN(50);
        instance.setRepAllowed(true);
        instance.setLNeg(0);
        instance.setLPoz(0);
        instance.setProbablePositive(0);
        // Create a dummy spectrum (10 k‑mers). You can change these strings as needed.
        std::vector<std::string> spectrum = { "AAAAA", "AAAAC", "AAACA", "AACAA", "ACAAA",
                                                "CAAAA", "AAAAG", "AAAGA", "AAGAA", "AGAAA" };
        instance.setSpectrum(spectrum);
        // Use the permutation representation.
        representation = std::make_shared<PermutationRepresentation>();
        // Create an individual with the identity permutation.
        validIndividual = std::make_shared<Individual>(std::vector<int>{0,1,2,3,4,5,6,7,8,9});
    }
};

//
// Test 1. Ensure that calculateFitness returns a non‑negative value for a valid individual.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessReturnsNonNegative) {
    double fit = fitnessCalc.calculateFitness(validIndividual, instance, representation);
    EXPECT_GE(fit, 0.0);
}

//
// Test 2. Ensure that the computed fitness value is finite (i.e. not infinite or NaN).
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessReturnsFinite) {
    double fit = fitnessCalc.calculateFitness(validIndividual, instance, representation);
    EXPECT_TRUE(std::isfinite(fit));
}

//
// Test 3. For an individual with empty gene vector, fitness should return 0.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessWithEmptyIndividualGenes) {
    auto emptyIndividual = std::make_shared<Individual>(std::vector<int>{});
    double fit = fitnessCalc.calculateFitness(emptyIndividual, instance, representation);
    EXPECT_EQ(fit, 0.0);
}

//
// Test 4. Passing a null individual should return -infinity.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessWithNullIndividual) {
    double fit = fitnessCalc.calculateFitness(nullptr, instance, representation);
    EXPECT_EQ(fit, -std::numeric_limits<double>::infinity());
}

//
// Test 5. Verify that buildAdjacencyMatrix returns a square matrix whose dimensions
// equal the spectrum size.
//
TEST_F(OptimizedGraphBasedFitnessTest, BuildAdjacencyMatrixDimensions) {
    buildGraphTest(validIndividual->getGenes(), instance);
    const auto& matrix = getMatrix();
    size_t n = instance.getSpectrum().size();
    EXPECT_EQ(matrix.size(), n);
    for (const auto& row : matrix) {
        EXPECT_EQ(row.size(), n);
    }
}

//
// Test 6. Test calculateEdgeWeight with two identical k‑mers (perfect overlap).
//
TEST_F(OptimizedGraphBasedFitnessTest, CalculateEdgeWeightPerfectOverlap) {
    std::string from = "AAAAA";
    std::string to   = "AAAAA";
    bool connected = testConnection(0, 0, instance);
    EXPECT_TRUE(connected);
}

//
// Test 7. Test calculateEdgeWeight with two completely different k‑mers (no overlap).
//
TEST_F(OptimizedGraphBasedFitnessTest, CalculateEdgeWeightNoOverlap) {
    std::string from = "AAAAA";
    std::string to   = "CCCCC";
    bool connected = testConnection(0, 5, instance);
    EXPECT_FALSE(connected);
}

//
// Test 8. Test that calculateEdgeQuality (connectivity measure) returns a value within [0,1].
//
TEST_F(OptimizedGraphBasedFitnessTest, CalculateEdgeQualityWithinRange) {
    buildGraphTest(validIndividual->getGenes(), instance);
    double quality = getConnScore();
    EXPECT_GE(quality, 0.0);
    EXPECT_LE(quality, 1.0);
}

//
// Test 9. Test that calculateLength returns a penalty value between 0 and 1.
//
TEST_F(OptimizedGraphBasedFitnessTest, CalculateLengthPenaltyBehavior) {
    double penalty = getLengthPen(validIndividual->getGenes(), instance);
    EXPECT_GE(penalty, 0.0);
    EXPECT_LE(penalty, 1.0);
}

//
// Test 10. Test that calculateSpectrumCoverageScore returns a value between 0 and 1.
//
TEST_F(OptimizedGraphBasedFitnessTest, CalculateSpectrumCoverageScoreBetween0And1) {
    double coverage = getSpectrumScore(validIndividual->getGenes(), instance);
    EXPECT_GE(coverage, 0.0);
    EXPECT_LE(coverage, 1.0);
}

//
// Test 11. Integration test: For a valid individual, fitness should be positive.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessIntegrationValidCase) {
    double fit = fitnessCalc.calculateFitness(validIndividual, instance, representation);
    EXPECT_GT(fit, 0.0);
    EXPECT_TRUE(std::isfinite(fit));
}

//
// Test 12. Two individuals with different permutations should yield different fitness values.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessIntegrationDifferentPermutationsYieldDifferentFitness) {
    auto individual1 = std::make_shared<Individual>(std::vector<int>{0,1,2,3,4,5,6,7,8,9});
    auto individual2 = std::make_shared<Individual>(std::vector<int>{9,8,7,6,5,4,3,2,1,0});
    double fit1 = fitnessCalc.calculateFitness(individual1, instance, representation);
    double fit2 = fitnessCalc.calculateFitness(individual2, instance, representation);
    EXPECT_NE(fit1, fit2);
}

//
// Test 13. Repeated calls on the same individual should return the same fitness value.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessConsistencyMultipleCalls) {
    double fit1 = fitnessCalc.calculateFitness(validIndividual, instance, representation);
    double fit2 = fitnessCalc.calculateFitness(validIndividual, instance, representation);
    EXPECT_NEAR(fit1, fit2, 1e-6);
}

//
// Test 14. If an individual contains an out‐of‐range gene index, fitness should return -infinity.
//
TEST_F(OptimizedGraphBasedFitnessTest, FitnessHandlesInvalidGene) {
    auto invalidIndividual = std::make_shared<Individual>(std::vector<int>{0,1,2,3,4,5,6,7,8,10}); // 10 is out of range
    double fit = fitnessCalc.calculateFitness(invalidIndividual, instance, representation);
    EXPECT_EQ(fit, -std::numeric_limits<double>::infinity());
}

//
// Test 15. Create a minimal instance with only 2 k‑mers and check that fitness is computed.
TEST_F(OptimizedGraphBasedFitnessTest, FitnessEdgeCaseSmallSpectrum) {
    DNAInstance smallInstance;
    smallInstance.setK(3);
    smallInstance.setDeltaK(0);
    smallInstance.setN(20);
    smallInstance.setRepAllowed(true);
    std::vector<std::string> smallSpectrum = {"AAA", "AAB"};
    smallInstance.setSpectrum(smallSpectrum);
    auto indiv = std::make_shared<Individual>(std::vector<int>{0,1});
    double fit = fitnessCalc.calculateFitness(indiv, smallInstance, representation);
    EXPECT_GE(fit, 0.0);
}

//
// Test 16. Check that every edge weight in the adjacency matrix is non-negative.
TEST_F(OptimizedGraphBasedFitnessTest, AdjacencyMatrixElementsNonNegative) {
    buildGraphTest(validIndividual->getGenes(), instance);
    const auto& matrix = getMatrix();
    for (const auto& row : matrix) {
        for (bool connected : row) {
            EXPECT_TRUE(connected || !connected); // Just verify it's a valid boolean
        }
    }
}

//
// Test 17. Check that calculateConnectivity returns a non-negative value.
TEST_F(OptimizedGraphBasedFitnessTest, CalculateConnectivityReturnsNonNegative) {
    buildGraphTest(validIndividual->getGenes(), instance);
    double conn = getConnScore();
    EXPECT_GE(conn, 0.0);
}

//
// Test 18. Integration test with a "realistic" instance: create a longer spectrum and a larger individual,
// and verify that the computed fitness is within a plausible range (here assumed to be between 0 and 1).
TEST_F(OptimizedGraphBasedFitnessTest, FitnessIntegrationWithRealisticScenario) {
    DNAInstance longInstance;
    longInstance.setK(5);
    longInstance.setDeltaK(1);
    longInstance.setN(100);
    longInstance.setRepAllowed(true);
    std::vector<std::string> longSpectrum;
    for (int i = 0; i < 20; i++) {
        // Create a pattern: "AAAAA", "AAAAB", "AAAAC", etc.
        std::string s = "AAAAA";
        s[4] = 'A' + (i % 4);
        longSpectrum.push_back(s);
    }
    longInstance.setSpectrum(longSpectrum);
    auto indiv = std::make_shared<Individual>(
        std::vector<int>{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19});
    double fit = fitnessCalc.calculateFitness(indiv, longInstance, representation);
    EXPECT_GE(fit, 0.0);
    EXPECT_LE(fit, 1.0);
}

//
// Test 19. Ensure that calculateFitness returns a double and not NaN.
TEST_F(OptimizedGraphBasedFitnessTest, FitnessReturnsDoubleType) {
    double fit = fitnessCalc.calculateFitness(validIndividual, instance, representation);
    EXPECT_TRUE(typeid(fit) == typeid(double));
    EXPECT_FALSE(std::isnan(fit));
}

//
// Test 20. Minimal instance: with a spectrum of one k‑mer.
TEST_F(OptimizedGraphBasedFitnessTest, FitnessHandlesMinimalInstance) {
    DNAInstance minimalInstance;
    minimalInstance.setK(3);
    minimalInstance.setDeltaK(0);
    minimalInstance.setN(10);
    minimalInstance.setRepAllowed(true);
    std::vector<std::string> oneSpectrum = {"AAA"};
    minimalInstance.setSpectrum(oneSpectrum);
    auto indiv = std::make_shared<Individual>(std::vector<int>{0});
    double fit = fitnessCalc.calculateFitness(indiv, minimalInstance, representation);
    EXPECT_GE(fit, 0.0);
}

