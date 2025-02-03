//
// Created by konrad_guest on 03/02/2025.
//
// OrderCrossover_test.cpp
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

// Helper: create a dummy DNA instance with a fixed spectrum.
DNAInstance createDummyDNAInstance() {
    DNAInstance instance;
    instance.setK(5);
    // For testing purposes, we create a spectrum of 10 fixed "k-mers"
    std::vector<std::string> spectrum = { "AAAAA", "AAAAC", "AAACA", "AACAA", "ACAAA",
                                            "CAAAA", "AAAAG", "AAAGA", "AAGAA", "AGAAA" };
    instance.setSpectrum(spectrum);
    instance.setN(50);
    instance.setDeltaK(0);
    instance.setRepAllowed(true);
    return instance;
}

// Helper: create an Individual from a given permutation vector.
std::shared_ptr<Individual> createIndividual(const std::vector<int>& perm) {
    return std::make_shared<Individual>(perm);
}

// Test fixture for OrderCrossover.
class OrderCrossoverTest : public BaseTest {
protected:
    DNAInstance instance;
    std::shared_ptr<IRepresentation> representation;
    static bool s_loggerInitialized;

    static void SetUpTestSuite() {
        if (!s_loggerInitialized) {
            try {
                Logger::initialize("order_crossover_test.log");
                Logger::setLogLevel(LogLevel::DEBUG);
                s_loggerInitialized = true;
                LOG_INFO("OrderCrossoverTest test suite initialized successfully");
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize logger: " << e.what() << std::endl;
                throw;
            }
        }
    }

    static void TearDownTestSuite() {
        if (s_loggerInitialized) {
            LOG_INFO("Tearing down OrderCrossoverTest test suite");
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
            
            LOG_DEBUG("OrderCrossoverTest::SetUp - Starting setup");
            BaseTest::SetUp();  // Call parent's SetUp first
            
            instance = createDummyDNAInstance();
            representation = std::make_shared<PermutationRepresentation>();
            
            LOG_DEBUG("OrderCrossoverTest::SetUp - Setup completed successfully");
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
            LOG_DEBUG("OrderCrossoverTest::TearDown - Starting cleanup");
            BaseTest::TearDown();  // Call parent's TearDown
            LOG_DEBUG("OrderCrossoverTest::TearDown - Cleanup completed successfully");
        } catch (const std::exception& e) {
            std::cerr << "TearDown failed: " << e.what() << std::endl;
            LOG_ERROR("TearDown failed: {}", e.what());
            throw;
        }
    }
};

bool OrderCrossoverTest::s_loggerInitialized = false;

TEST_F(OrderCrossoverTest, OffspringCountNonZero) {
    OrderCrossover oc;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividual({0,1,2,3,4,5,6,7,8,9}));
    parents.push_back(createIndividual({9,8,7,6,5,4,3,2,1,0}));
    auto offspring = oc.crossover(parents, instance, representation);
    EXPECT_FALSE(offspring.empty());
    EXPECT_GE(offspring.size(), 1);
    EXPECT_LE(offspring.size(), 2);
}

TEST_F(OrderCrossoverTest, OffspringAreValidPermutation) {
    OrderCrossover oc;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividual({0,1,2,3,4,5,6,7,8,9}));
    parents.push_back(createIndividual({9,8,7,6,5,4,3,2,1,0}));
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        const auto& genes = child->getGenes();
        std::vector<bool> seen(10, false);
        for (int gene : genes) {
            EXPECT_GE(gene, 0);
            EXPECT_LT(gene, 10);
            EXPECT_FALSE(seen[gene]);
            seen[gene] = true;
        }
    }
}

TEST_F(OrderCrossoverTest, OffspringGeneSetSameAsParents) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = {createIndividual(parentGenes1), createIndividual(parentGenes2)};
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        auto genes = child->getGenes();
        std::sort(genes.begin(), genes.end());
        std::vector<int> expected = {0,1,2,3,4,5,6,7,8,9};
        EXPECT_EQ(genes, expected);
    }
}

TEST_F(OrderCrossoverTest, RepeatedCrossoverProducesDifferentOffspring) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes1), createIndividual(parentGenes2) };
    std::set<std::vector<int>> offspringSet;
    for (int i = 0; i < 50; ++i) {
        auto offspring = oc.crossover(parents, instance, representation);
        for (auto& child : offspring) {
            offspringSet.insert(child->getGenes());
        }
    }
    EXPECT_GT(offspringSet.size(), 1);
}

TEST_F(OrderCrossoverTest, OffspringNotEqualToAnyParent) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {0,2,4,6,8,9,7,5,3,1};
    auto parent1 = createIndividual(parentGenes1);
    auto parent2 = createIndividual(parentGenes2);
    std::vector<std::shared_ptr<Individual>> parents = {parent1, parent2};
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        bool equalToParent1 = (child->getGenes() == parentGenes1);
        bool equalToParent2 = (child->getGenes() == parentGenes2);
        EXPECT_FALSE(equalToParent1 && equalToParent2);
    }
}

TEST_F(OrderCrossoverTest, PreservesSegmentOrder) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {3,4,5,6,7,8,9,0,1,2};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes1), createIndividual(parentGenes2) };
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        auto genes = child->getGenes();
        bool found = false;
        for (size_t i = 0; i + 2 < genes.size(); ++i) {
            std::vector<int> sub(genes.begin() + i, genes.begin() + i + 3);
            if (std::search(parentGenes1.begin(), parentGenes1.end(), sub.begin(), sub.end()) != parentGenes1.end()) {
                found = true;
                break;
            }
        }
        EXPECT_TRUE(found);
    }
}

TEST_F(OrderCrossoverTest, SingleParentEdgeCase) {
    OrderCrossover oc;
    std::vector<std::shared_ptr<Individual>> parents;
    parents.push_back(createIndividual({0,1,2,3,4,5,6,7,8,9}));
    auto offspring = oc.crossover(parents, instance, representation);
    EXPECT_EQ(offspring.size(), 1);
    EXPECT_EQ(offspring[0]->getGenes(), (std::vector<int>{0,1,2,3,4,5,6,7,8,9}));
}

TEST_F(OrderCrossoverTest, EmptyParentsVector) {
    OrderCrossover oc;
    std::vector<std::shared_ptr<Individual>> parents;
    auto offspring = oc.crossover(parents, instance, representation);
    EXPECT_TRUE(offspring.empty());
}

TEST_F(OrderCrossoverTest, InvalidParentGenes) {
    OrderCrossover oc;
    std::vector<std::shared_ptr<Individual>> parents;
    // Introduce an invalid gene (10 is out of range for our spectrum of size 10)
    parents.push_back(createIndividual({0,1,2,3,4,10,6,7,8,9}));
    parents.push_back(createIndividual({9,8,7,6,5,4,3,2,1,0}));
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_FALSE(representation->isValid(child, instance));
    }
}

TEST_F(OrderCrossoverTest, OffspringHaveCorrectSize) {
    OrderCrossover oc;
    std::vector<int> parentGenes = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes), createIndividual(parentGenes) };
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        EXPECT_EQ(child->getGenes().size(), parentGenes.size());
    }
}

TEST_F(OrderCrossoverTest, MultipleCallsConsistency) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes1), createIndividual(parentGenes2) };
    auto offspring1 = oc.crossover(parents, instance, representation);
    auto offspring2 = oc.crossover(parents, instance, representation);
    ASSERT_EQ(offspring1.size(), offspring2.size());
    for (size_t i = 0; i < offspring1.size(); i++) {
        EXPECT_EQ(offspring1[i]->getGenes().size(), offspring2[i]->getGenes().size());
    }
}

TEST_F(OrderCrossoverTest, OffspringFitnessIsComputable) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes1), createIndividual(parentGenes2) };
    auto offspring = oc.crossover(parents, instance, representation);
    OptimizedGraphBasedFitness fitnessCalc;
    for (auto& child : offspring) {
        double fit = fitnessCalc.calculateFitness(child, instance, representation);
        EXPECT_GE(fit, 0.0);
    }
}

TEST_F(OrderCrossoverTest, OffspringGeneOrderRandomness) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {1,2,3,4,5,6,7,8,9,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes1), createIndividual(parentGenes2) };
    std::vector<std::vector<int>> offspringGenes;
    for (int i = 0; i < 30; i++) {
        auto offspring = oc.crossover(parents, instance, representation);
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

TEST_F(OrderCrossoverTest, DNAAssemblyFromOffspring) {
    OrderCrossover oc;
    std::vector<int> parentGenes1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> parentGenes2 = {9,8,7,6,5,4,3,2,1,0};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes1), createIndividual(parentGenes2) };
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        std::string dna = representation->toDNA(child, instance);
        EXPECT_GE(dna.length(), instance.getK());
        EXPECT_LE(dna.length(), instance.getN() * 2);
    }
}

TEST_F(OrderCrossoverTest, OffspringStabilityUnderRepeatedCrossover) {
    OrderCrossover oc;
    std::vector<int> parentGenes = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes), createIndividual(parentGenes) };
    auto offspring1 = oc.crossover(parents, instance, representation);
    auto offspring2 = oc.crossover(parents, instance, representation);
    for (auto& child : offspring1) {
        EXPECT_EQ(child->getGenes(), parentGenes);
    }
    for (auto& child : offspring2) {
        EXPECT_EQ(child->getGenes(), parentGenes);
    }
}

TEST_F(OrderCrossoverTest, TestOffspringInputIntegrity) {
    OrderCrossover oc;
    std::vector<int> perm1 = {0,1,2,3,4,5,6,7,8,9};
    std::vector<int> perm2 = {9,8,7,6,5,4,3,2,1,0};
    auto parent1 = createIndividual(perm1);
    auto parent2 = createIndividual(perm2);
    std::vector<std::shared_ptr<Individual>> parents = { parent1, parent2 };
    oc.crossover(parents, instance, representation);
    EXPECT_EQ(parent1->getGenes(), perm1);
    EXPECT_EQ(parent2->getGenes(), perm2);
}

TEST_F(OrderCrossoverTest, TestOffspringGeneRangeWithinSpectrum) {
    OrderCrossover oc;
    std::vector<int> parentGenes = {0,1,2,3,4,5,6,7,8,9};
    std::vector<std::shared_ptr<Individual>> parents = { createIndividual(parentGenes), createIndividual(parentGenes) };
    auto offspring = oc.crossover(parents, instance, representation);
    for (auto& child : offspring) {
        for (int gene : child->getGenes()) {
            EXPECT_GE(gene, 0);
            EXPECT_LT(gene, static_cast<int>(instance.getSpectrum().size()));
        }
    }
}
