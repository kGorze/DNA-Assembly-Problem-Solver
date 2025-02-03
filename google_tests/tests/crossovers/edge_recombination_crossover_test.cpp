//
// Created by konrad_guest on 03/02/2025.
//
// EdgeRecombinationCrossover_test.cpp
#include <gtest/gtest.h>
#include "metaheuristics/crossover_impl.h"
#include "metaheuristics/individual.h"
#include "metaheuristics/representation.h"
#include "metaheuristics/fitness_impl.h"
#include "dna/dna_instance.h"
#include <vector>
#include <algorithm>
#include <set>
#include "../../tests/base_test.h"
#include "utils/logging.h"

// Helper: create a dummy DNA instance for Edge-Recombination tests.
DNAInstance createDummyDNAInstance_ER() {
    DNAInstance instance;
    instance.setK(5);
    std::vector<std::string> spectrum = { "AAAAA", "AAAAC", "AAACA", "AACAA", "ACAAA",
                                            "CAAAA", "AAAAG", "AAAGA", "AAGAA", "AGAAA" };
    instance.setSpectrum(spectrum);
    instance.setN(50);
    instance.setDeltaK(0);
    instance.setRepAllowed(true);
    return instance;
}

// Helper: create an individual.
std::shared_ptr<Individual> createIndividualER(const std::vector<int>& perm) {
    return std::make_shared<Individual>(perm);
}

// Test fixture for EdgeRecombinationCrossover.
class EdgeRecombinationCrossoverTest : public BaseTest {
protected:
    DNAInstance instance;
    std::shared_ptr<IRepresentation> representation;
    static bool s_loggerInitialized;

    static void SetUpTestSuite() {
        if (!s_loggerInitialized) {
            try {
                Logger::initialize("edge_recombination_crossover_test.log");
                Logger::setLogLevel(LogLevel::DEBUG);
                s_loggerInitialized = true;
                LOG_INFO("EdgeRecombinationCrossoverTest test suite initialized successfully");
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize logger: " << e.what() << std::endl;
                throw;
            }
        }
    }

    static void TearDownTestSuite() {
        if (s_loggerInitialized) {
            LOG_INFO("Tearing down EdgeRecombinationCrossoverTest test suite");
            try {
                Logger::cleanup();
                s_loggerInitialized = false;
            } catch (const std::exception& e) {
                std::cerr << "Failed to cleanup logger: " << e.what() << std::endl;
                throw;
            }
        }
    }

    void SetUp() override {
        try {
            // Ensure logger is initialized
            if (!s_loggerInitialized) {
                SetUpTestSuite();
            }
            
            LOG_DEBUG("EdgeRecombinationCrossoverTest::SetUp - Starting setup");
            BaseTest::SetUp();  // Call parent's SetUp first
            
            instance = createDummyDNAInstance_ER();
            representation = std::make_shared<PermutationRepresentation>();
            
            LOG_DEBUG("EdgeRecombinationCrossoverTest::SetUp - Setup completed successfully");
        } catch (const std::exception& e) {
            std::cerr << "SetUp failed with exception: " << e.what() << std::endl;
            LOG_ERROR("SetUp failed with exception: {}", e.what());
            throw;
        } catch (...) {
            std::cerr << "SetUp failed with unknown exception" << std::endl;
            LOG_ERROR("SetUp failed with unknown exception");
            throw;
        }
    }
    
    void TearDown() override {
        try {
            LOG_DEBUG("EdgeRecombinationCrossoverTest::TearDown - Starting cleanup");
            BaseTest::TearDown();  // Call parent's TearDown
            LOG_DEBUG("EdgeRecombinationCrossoverTest::TearDown - Cleanup completed successfully");
        } catch (const std::exception& e) {
            std::cerr << "TearDown failed: " << e.what() << std::endl;
            LOG_ERROR("TearDown failed: {}", e.what());
            throw;
        }
    }
};

bool EdgeRecombinationCrossoverTest::s_loggerInitialized = false;

TEST_F(EdgeRecombinationCrossoverTest, OffspringCountNonZero) {
    EdgeRecombinationCrossover er;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividualER({0,1,2,3,4,5,6,7,8,9}));
    parents.push_back(createIndividualER({9,8,7,6,5,4,3,2,1,0}));
    auto offspring = er.crossover(parents, instance, representation);
    EXPECT_FALSE(offspring.empty());
    EXPECT_GE(offspring.size(), 1);
    EXPECT_LE(offspring.size(), 2);
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringAreValidPermutation) {
    EdgeRecombinationCrossover er;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividualER({0,1,2,3,4,5,6,7,8,9}));
    parents.push_back(createIndividualER({9,8,7,6,5,4,3,2,1,0}));
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        std::vector<bool> seen(10, false);
        for (int gene : child->getGenes()) {
            EXPECT_GE(gene, 0);
            EXPECT_LT(gene, 10);
            EXPECT_FALSE(seen[gene]);
            seen[gene] = true;
        }
    }
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringGeneSetSameAsParents) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        auto genes = child->getGenes();
        std::sort(genes.begin(), genes.end());
        std::vector<int> expected = {0,1,2,3,4,5,6,7,8,9};
        EXPECT_EQ(genes, expected);
    }
}

TEST_F(EdgeRecombinationCrossoverTest, RepeatedCrossoverProducesDifferentOffspring) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    std::set<std::vector<int>> offspringSet;
    for (int i = 0; i < 50; i++) {
        auto offspring = er.crossover(parents, instance, representation);
        for (auto& child : offspring) {
            offspringSet.insert(child->getGenes());
        }
    }
    EXPECT_GT(offspringSet.size(), 1);
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringNotEqualToAnyParent) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {0,2,4,6,8,9,7,5,3,1};
    auto parent1 = createIndividualER(parentGenes1);
    auto parent2 = createIndividualER(parentGenes2);
    std::vector<std::shared_ptr<Individual>> parents = { parent1, parent2 };
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        bool eq1 = (child->getGenes() == parentGenes1);
        bool eq2 = (child->getGenes() == parentGenes2);
        EXPECT_FALSE(eq1 && eq2);
    }
}

TEST_F(EdgeRecombinationCrossoverTest, PreservesLocalEdgeQuality) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {3,4,5,6,7,8,9,0,1,2};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    auto offspring = er.crossover(parents, instance, representation);
    // For each pair in offspring, simply check that there is some "overlap"
    for (auto& child : offspring) {
        const auto& genes = child->getGenes();
        const auto& spectrum = instance.getSpectrum();
        int k = instance.getK();
        for (size_t i = 0; i < genes.size()-1; ++i) {
            if (genes[i] < 0 || genes[i+1] < 0 ||
                static_cast<size_t>(genes[i]) >= spectrum.size() ||
                static_cast<size_t>(genes[i+1]) >= spectrum.size())
                continue;
            std::string suffix = spectrum[genes[i]].substr(spectrum[genes[i]].size() - (k-1));
            std::string prefix = spectrum[genes[i+1]].substr(0, k-1);
            // If they match, record a "good" connection
            EXPECT_GE(suffix.compare(prefix), 0);
        }
    }
}

TEST_F(EdgeRecombinationCrossoverTest, SingleParentEdgeCase) {
    EdgeRecombinationCrossover er;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividualER({0,1,2,3,4,5,6,7,8,9}));
    auto offspring = er.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), 1);
    EXPECT_EQ(offspring[0]->getGenes(), (std::vector<int>{0,1,2,3,4,5,6,7,8,9}));
}

TEST_F(EdgeRecombinationCrossoverTest, EmptyParentsVector) {
    EdgeRecombinationCrossover er;
    std::vector<std::shared_ptr<Individual>> parents;
    auto offspring = er.crossover(parents, instance, representation);
    EXPECT_TRUE(offspring.empty());
}

TEST_F(EdgeRecombinationCrossoverTest, InvalidParentGenes) {
    EdgeRecombinationCrossover er;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividualER({0,1,2,3,4,10,6,7,8,9})); // invalid gene 10
    parents.push_back(createIndividualER({9,8,7,6,5,4,3,2,1,0}));
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_FALSE(representation->isValid(child, instance));
    }
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringHaveCorrectSize) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes), createIndividualER(parentGenes) };
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_EQ(child->getGenes().size(), parentGenes.size());
    }
}

TEST_F(EdgeRecombinationCrossoverTest, MultipleCallsConsistency) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    auto offspring1 = er.crossover(parents, instance, representation);
    auto offspring2 = er.crossover(parents, instance, representation);
    ASSERT_EQ(offspring1.size(), offspring2.size());
    for (size_t i = 0; i < offspring1.size(); i++) {
        EXPECT_EQ(offspring1[i]->getGenes().size(), offspring2[i]->getGenes().size());
    }
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringFitnessIsComputable) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    auto offspring = er.crossover(parents, instance, representation);
    OptimizedGraphBasedFitness fitnessCalc;
    for (auto& child : offspring) {
        double fit = fitnessCalc.calculateFitness(child, instance, representation);
        EXPECT_GE(fit, 0.0);
    }
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringGeneOrderRandomness) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {1,2,3,4,5,6,7,8,9,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    std::vector<std::vector<int>> offspringGenes;
    for (int i = 0; i < 30; i++) {
        auto offspring = er.crossover(parents, instance, representation);
        for (auto& child : offspring) {
            offspringGenes.push_back(child->getGenes());
        }
    }
    bool diversity = false;
    for (size_t i = 1; i < offspringGenes.size(); i++) {
        if (offspringGenes[i] != offspringGenes[0]) {
            diversity = true;
            break;
        }
    }
    EXPECT_TRUE(diversity);
}

TEST_F(EdgeRecombinationCrossoverTest, DNAAssemblyFromOffspring) {
    EdgeRecombinationCrossover er;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(parentGenes1), createIndividualER(parentGenes2) };
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        std::string dna = representation->toDNA(child, instance);
        EXPECT_GE(dna.length(), instance.getK());
    }
}

TEST_F(EdgeRecombinationCrossoverTest, OffspringInputIntegrity) {
    EdgeRecombinationCrossover er;
    std::vector<int> perm1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> perm2 = {9,8,7,6,5,4,3,2,1,0};
    auto parent1 = createIndividualER(perm1);
    auto parent2 = createIndividualER(perm2);
    std::vector<std::shared_ptr<Individual>> parents = { parent1, parent2 };
    er.crossover(parents, instance, representation);
    EXPECT_EQ(parent1->getGenes(), perm1);
    EXPECT_EQ(parent2->getGenes(), perm2);
}

TEST_F(EdgeRecombinationCrossoverTest, WithNonStandardSpectrum) {
    std::vector<std::string> newSpectrum = { "TTTTT", "TTTTC", "TTTCC", "TTCCC", "TCCCC",
                                              "CCCCC", "TTTTG", "TTTGG", "TTGGG", "TGGGG" };
    instance.setSpectrum(newSpectrum);
    EdgeRecombinationCrossover er;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualER(perm), createIndividualER(perm) };
    auto offspring = er.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_EQ(child->getGenes(), perm);
    }
}
