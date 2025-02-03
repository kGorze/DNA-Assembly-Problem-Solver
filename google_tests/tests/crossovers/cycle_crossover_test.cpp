//
// Created by konrad_guest on 03/02/2025.
//
// CycleCrossover_test.cpp
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

// Helper: create a dummy DNA instance for CycleCrossover tests.
DNAInstance createDummyDNAInstance_Cycle() {
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
std::shared_ptr<Individual> createIndividualCycle(const std::vector<int>& perm) {
    return std::make_shared<Individual>(perm);
}

// Test fixture for CycleCrossover.
class CycleCrossoverTest : public BaseTest {
protected:
    DNAInstance instance;
    std::shared_ptr<IRepresentation> representation;
    static bool s_loggerInitialized;

    static void SetUpTestSuite() {
        if (!s_loggerInitialized) {
            try {
                Logger::initialize("cycle_crossover_test.log");
                Logger::setLogLevel(LogLevel::DEBUG);
                s_loggerInitialized = true;
                LOG_INFO("CycleCrossoverTest test suite initialized successfully");
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize logger: " << e.what() << std::endl;
                throw;
            }
        }
    }

    static void TearDownTestSuite() {
        if (s_loggerInitialized) {
            LOG_INFO("Tearing down CycleCrossoverTest test suite");
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
            
            LOG_DEBUG("CycleCrossoverTest::SetUp - Starting setup");
            BaseTest::SetUp();  // Call parent's SetUp first
            
            instance = createDummyDNAInstance_Cycle();
            representation = std::make_shared<PermutationRepresentation>();
            
            LOG_DEBUG("CycleCrossoverTest::SetUp - Setup completed successfully");
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
            LOG_DEBUG("CycleCrossoverTest::TearDown - Starting cleanup");
            BaseTest::TearDown();  // Call parent's TearDown
            LOG_DEBUG("CycleCrossoverTest::TearDown - Cleanup completed successfully");
        } catch (const std::exception& e) {
            std::cerr << "TearDown failed: " << e.what() << std::endl;
            LOG_ERROR("TearDown failed: {}", e.what());
            throw;
        }
    }
};

bool CycleCrossoverTest::s_loggerInitialized = false;

TEST_F(CycleCrossoverTest, ReturnsParentsAsOffspring) {
    CycleCrossover cc;
    std::vector<std::shared_ptr<Individual>> parents;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    parents.push_back(createIndividualCycle(perm));
    std::vector<int> rev = perm; std::reverse(rev.begin(), rev.end());
    parents.push_back(createIndividualCycle(rev));
    auto offspring = cc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), parents.size());
    for (size_t i = 0; i < parents.size(); i++) {
        EXPECT_EQ(offspring[i]->getGenes(), parents[i]->getGenes());
    }
}

TEST_F(CycleCrossoverTest, HandlesEmptyParents) {
    CycleCrossover cc;
    std::vector<std::shared_ptr<Individual>> parents;
    auto offspring = cc.crossover(parents, instance, representation);
    EXPECT_TRUE(offspring.empty());
}

TEST_F(CycleCrossoverTest, HandlesSingleParent) {
    CycleCrossover cc;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividualCycle({0,1,2,3,4,5,6,7,8,9}));
    auto offspring = cc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), 1);
    EXPECT_EQ(offspring[0]->getGenes(), (std::vector<int>{0,1,2,3,4,5,6,7,8,9}));
}

TEST_F(CycleCrossoverTest, OffspringAreValidPermutation) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, representation);
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

TEST_F(CycleCrossoverTest, OffspringGeneSetSameAsParents) {
    CycleCrossover cc;
    std::vector<int> perm1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> perm2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm1), createIndividualCycle(perm2) };
    auto offspring = cc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), 2);
    EXPECT_EQ(offspring[0]->getGenes(), perm1);
    EXPECT_EQ(offspring[1]->getGenes(), perm2);
}

TEST_F(CycleCrossoverTest, RepeatedCallsReturnSameParents) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring1 = cc.crossover(parents, instance, representation);
    auto offspring2 = cc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring1.size(), offspring2.size());
    for (size_t i = 0; i < offspring1.size(); i++) {
        EXPECT_EQ(offspring1[i]->getGenes(), offspring2[i]->getGenes());
    }
}

TEST_F(CycleCrossoverTest, OffspringDNAAssemblyMatchesParent) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        std::string dna = representation->toDNA(child, instance);
        std::string expected = representation->toDNA(parents[0], instance);
        EXPECT_EQ(dna, expected);
    }
}

TEST_F(CycleCrossoverTest, TestWithDifferentParentsButCycleReturnsParents) {
    CycleCrossover cc;
    std::vector<int> perm1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> perm2 = {1,2,3,4,5,6,7,8,9,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm1), createIndividualCycle(perm2) };
    auto offspring = cc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), 2);
    EXPECT_EQ(offspring[0]->getGenes(), perm1);
    EXPECT_EQ(offspring[1]->getGenes(), perm2);
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverWithIdenticalParents) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_EQ(child->getGenes(), perm);
    }
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverPerformanceNotCrashing) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    for (int i = 0; i < 50; i++) {
        auto offspring = cc.crossover(parents, instance, representation);
        EXPECT_EQ(offspring.size(), 2);
    }
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverHandlesLargePermutation) {
    CycleCrossover cc;
    std::vector<int> largePerm(100);
    std::iota(largePerm.begin(), largePerm.end(), 0);
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(largePerm), createIndividualCycle(largePerm) };
    auto offspring = cc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), 2);
    for (auto& child : offspring) {
        EXPECT_EQ(child->getGenes().size(), largePerm.size());
    }
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverMultipleRandomRuns) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    std::set<std::vector<int>> uniqueResults;
    for (int i = 0; i < 30; i++) {
        auto offspring = cc.crossover(parents, instance, representation);
        for (auto& child : offspring) {
            uniqueResults.insert(child->getGenes());
        }
    }
    // Since cycle crossover simply returns the input, all offspring should be identical.
    EXPECT_EQ(uniqueResults.size(), 1);
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverDoesNotAlterInputParents) {
    CycleCrossover cc;
    std::vector<int> perm1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> perm2 = {9,8,7,6,5,4,3,2,1,0};
    auto parent1 = createIndividualCycle(perm1);
    auto parent2 = createIndividualCycle(perm2);
    std::vector<std::shared_ptr<Individual>> parents = { parent1, parent2 };
    cc.crossover(parents, instance, representation);
    EXPECT_EQ(parent1->getGenes(), perm1);
    EXPECT_EQ(parent2->getGenes(), perm2);
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverReturnSameForIdenticalParents) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_EQ(child->getGenes(), perm);
    }
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverInputIntegrity) {
    CycleCrossover cc;
    std::vector<int> perm1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> perm2 = {9,8,7,6,5,4,3,2,1,0};
    auto parent1 = createIndividualCycle(perm1);
    auto parent2 = createIndividualCycle(perm2);
    std::vector<std::shared_ptr<Individual>> parents = { parent1, parent2 };
    cc.crossover(parents, instance, representation);
    EXPECT_EQ(parent1->getGenes(), perm1);
    EXPECT_EQ(parent2->getGenes(), perm2);
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverGeneRangeWithinSpectrum) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        for (int gene : child->getGenes()) {
            EXPECT_GE(gene, 0);
            EXPECT_LT(gene, static_cast<int>(instance.getSpectrum().size()));
        }
    }
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverDNAAssembly) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        std::string dna = representation->toDNA(child, instance);
        EXPECT_GE(dna.length(), instance.getK());
    }
}

TEST_F(CycleCrossoverTest, TestCycleCrossoverHandlesNullRepresentationGracefully) {
    CycleCrossover cc;
    std::vector<int> perm = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(perm), createIndividualCycle(perm) };
    auto offspring = cc.crossover(parents, instance, nullptr);
    EXPECT_TRUE(offspring.empty());
}

TEST_F(CycleCrossoverTest, OffspringFitnessIsComputable) {
    CycleCrossover cc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividualCycle(parentGenes1), createIndividualCycle(parentGenes2) };
    auto offspring = cc.crossover(parents, instance, representation);
    OptimizedGraphBasedFitness fitnessCalc;
    for (auto& child : offspring) {
        double fit = fitnessCalc.calculateFitness(child, instance, representation);
        EXPECT_GE(fit, 0.0);
    }
}
