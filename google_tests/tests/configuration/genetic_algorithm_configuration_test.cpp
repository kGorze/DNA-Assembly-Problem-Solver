//
// Created by konrad_guest on 02/02/2025.
//
// GeneticAlgorithmConfiguration_test.cpp
#include <gtest/gtest.h>
#include <fstream>
#include <sstream>
#include <cstdio>    // for remove()
#include <set>
#include <algorithm>

// Include the configuration header and any dependent headers.
#include "configuration/genetic_algorithm_configuration.h"
#include "metaheuristics/population_cache_impl.h"
#include "tuning/parameters_parser.h"  // Assumed to define ParameterSet
#include "../base_test.h"
#include "utils/logging.h"

// Fixture for GA configuration tests.
class GeneticAlgorithmConfigurationTest : public BaseTest {
protected:
    // Each test gets its own configuration instance.
    GAConfig config;
    static bool s_loggerInitialized;

    static void SetUpTestSuite() {
        if (!s_loggerInitialized) {
            try {
                Logger::initialize("genetic_algorithm_configuration_test.log");
                Logger::setLogLevel(LogLevel::DEBUG);
                s_loggerInitialized = true;
                LOG_INFO("GeneticAlgorithmConfigurationTest test suite initialized successfully");
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize logger: " << e.what() << std::endl;
                throw;
            }
        }
    }

    static void TearDownTestSuite() {
        if (s_loggerInitialized) {
            LOG_INFO("Tearing down GeneticAlgorithmConfigurationTest test suite");
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
            
            LOG_DEBUG("GeneticAlgorithmConfigurationTest::SetUp - Starting setup");
            BaseTest::SetUp();  // Call parent's SetUp first
            
            // For tests that need a cache (like getFitness, getReplacement) we can set a dummy cache.
            // However, many tests will not call these functions.
            // You can choose to set it here or within the specific tests.
            // Here we leave it unset by default.
            LOG_DEBUG("GeneticAlgorithmConfigurationTest::SetUp - Setup completed successfully");
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
            LOG_DEBUG("GeneticAlgorithmConfigurationTest::TearDown - Starting cleanup");
            BaseTest::TearDown();  // Call parent's TearDown
            LOG_DEBUG("GeneticAlgorithmConfigurationTest::TearDown - Cleanup completed successfully");
        } catch (const std::exception& e) {
            std::cerr << "TearDown failed: " << e.what() << std::endl;
            LOG_ERROR("TearDown failed: {}", e.what());
            throw;
        }
    }

    // Helper function: create a temporary config file with given contents and return its filename.
    std::string createTempConfigFile(const std::string& contents) {
        std::string filename = "temp_test_config.cfg";
        std::ofstream out(filename);
        out << contents;
        out.close();
        return filename;
    }

    // Helper function: remove the temporary file.
    void removeTempFile(const std::string& filename) {
        std::remove(filename.c_str());
    }
};

bool GeneticAlgorithmConfigurationTest::s_loggerInitialized = false;

//
// 1. Test that the default configuration values match what we expect.
//
TEST_F(GeneticAlgorithmConfigurationTest, DefaultConfigurationValues) {
    EXPECT_EQ(config.getPopulationSize(), 100);
    EXPECT_NEAR(config.getMutationRate(), 0.1, 1e-6);
    EXPECT_NEAR(config.getCrossoverProbability(), 0.8, 1e-6);
    EXPECT_EQ(config.getTournamentSize(), 2);
    EXPECT_NEAR(config.getReplacementRatio(), 0.7, 1e-6);
    EXPECT_EQ(config.getNoImprovementGenerations(), 30);
    EXPECT_EQ(config.getTimeLimitSeconds(), 60);
    EXPECT_NEAR(config.getTargetFitness(), 1.0, 1e-6);
    EXPECT_EQ(config.getSelectionMethod(), "rank");
    // Instance-specific defaults (if not overwritten by config file):
    EXPECT_EQ(config.getK(), 7);
    EXPECT_EQ(config.getDeltaK(), 0);
    EXPECT_EQ(config.getLNeg(), 10);
    EXPECT_EQ(config.getLPoz(), 10);
    EXPECT_TRUE(config.isRepAllowed());
    EXPECT_EQ(config.getProbablePositive(), 0);
    // Adaptive parameters:
    EXPECT_TRUE(config.getAdaptiveParams().useAdaptiveMutation);
    EXPECT_NEAR(config.getAdaptiveParams().minMutationRate, 0.1, 1e-6);
    EXPECT_NEAR(config.getAdaptiveParams().maxMutationRate, 0.4, 1e-6);
    EXPECT_EQ(config.getAdaptiveParams().stagnationGenerations, 5);
    EXPECT_NEAR(config.getAdaptiveParams().improvementThreshold, 0.01, 1e-6);
    // Diversity parameters:
    EXPECT_NEAR(config.getDiversityParams().sharingRadius, 0.2, 1e-6);
    EXPECT_TRUE(config.getDiversityParams().useFitnessSharing);
    EXPECT_FALSE(config.getDiversityParams().useCrowding);
    EXPECT_NEAR(config.getDiversityParams().sharingAlpha, 1.0, 1e-6);
}

//
// 2. Test that resetToDefaults restores default values after changes.
//
TEST_F(GeneticAlgorithmConfigurationTest, ResetToDefaultsRestoresDefaults) {
    config.setPopulationSize(500);
    config.setMutationRate(0.9);
    config.setTournamentSize(10);
    config.setCrossoverProbability(0.5);
    config.resetToDefaults();
    EXPECT_EQ(config.getPopulationSize(), 100);
    EXPECT_NEAR(config.getMutationRate(), 0.1, 1e-6);
    EXPECT_EQ(config.getTournamentSize(), 2);
    EXPECT_NEAR(config.getCrossoverProbability(), 0.8, 1e-6);
}

//
// 3. Test loading a valid config file: we write a temporary file with a few key=value pairs and verify the config is updated.
//
TEST_F(GeneticAlgorithmConfigurationTest, LoadConfigFromValidFile) {
    std::string configText =
        "populationSize = 200\n"
        "mutationRate = 0.2\n"
        "crossoverProbability = 0.9\n"
        "tournamentSize = 5\n"
        "replacementRatio = 0.6\n"
        "noImprovementGenerations = 50\n"
        "timeLimitSeconds = 120\n"
        "targetFitness = 0.95\n"
        "selectionMethod = tournament\n"
        "k = 8\n"
        "deltaK = 1\n"
        "lNeg = 0\n"
        "lPoz = 0\n"
        "repAllowed = false\n"
        "probablePositive = 5\n"
        "useAdaptiveMutation = false\n"
        "minMutationRate = 0.15\n"
        "maxMutationRate = 0.35\n"
        "stagnationGenerations = 7\n"
        "improvementThreshold = 0.02\n"
        "useFitnessSharing = false\n"
        "useCrowding = true\n"
        "sharingRadius = 0.25\n"
        "sharingAlpha = 1.2\n"
        "diversityWeight = 0.5\n"
        "crossoverType = adaptive\n";
    std::string filename = createTempConfigFile(configText);

    bool loaded = config.loadFromFile(filename);
    EXPECT_TRUE(loaded);
    EXPECT_EQ(config.getPopulationSize(), 200);
    EXPECT_NEAR(config.getMutationRate(), 0.2, 1e-6);
    EXPECT_NEAR(config.getCrossoverProbability(), 0.9, 1e-6);
    EXPECT_EQ(config.getTournamentSize(), 5);
    EXPECT_NEAR(config.getReplacementRatio(), 0.6, 1e-6);
    EXPECT_EQ(config.getNoImprovementGenerations(), 50);
    EXPECT_EQ(config.getTimeLimitSeconds(), 120);
    EXPECT_NEAR(config.getTargetFitness(), 0.95, 1e-6);
    EXPECT_EQ(config.getSelectionMethod(), "tournament");
    EXPECT_EQ(config.getK(), 8);
    EXPECT_EQ(config.getDeltaK(), 1);
    EXPECT_EQ(config.getLNeg(), 0);
    EXPECT_EQ(config.getLPoz(), 0);
    EXPECT_FALSE(config.isRepAllowed());
    EXPECT_EQ(config.getProbablePositive(), 5);
    EXPECT_FALSE(config.getAdaptiveParams().useAdaptiveMutation);
    EXPECT_NEAR(config.getAdaptiveParams().minMutationRate, 0.15, 1e-6);
    EXPECT_NEAR(config.getAdaptiveParams().maxMutationRate, 0.35, 1e-6);
    EXPECT_EQ(config.getAdaptiveParams().stagnationGenerations, 7);
    EXPECT_NEAR(config.getAdaptiveParams().improvementThreshold, 0.02, 1e-6);
    EXPECT_FALSE(config.getDiversityParams().useFitnessSharing);
    EXPECT_TRUE(config.getDiversityParams().useCrowding);
    EXPECT_NEAR(config.getDiversityParams().sharingRadius, 0.25, 1e-6);
    EXPECT_NEAR(config.getDiversityParams().sharingAlpha, 1.2, 1e-6);
    EXPECT_NEAR(config.getDiversityParams().diversityWeight, 0.5, 1e-6);
    // Clean up temporary file.
    removeTempFile(filename);
}

//
// 4. Test that loading a non-existent config file returns false.
//
TEST_F(GeneticAlgorithmConfigurationTest, LoadNonexistentConfigFileReturnsFalse) {
    bool loaded = config.loadFromFile("this_file_should_not_exist.cfg");
    EXPECT_FALSE(loaded);
}

//
// 5. Test that an invalid population size causes validate() to fail.
//
TEST_F(GeneticAlgorithmConfigurationTest, ValidateFailsWithInvalidPopulationSize) {
    config.setPopulationSize(0);
    EXPECT_FALSE(config.validate());
}

//
// 6. Test that an invalid mutation rate (e.g. 1.5) causes validate() to fail.
//
TEST_F(GeneticAlgorithmConfigurationTest, ValidateFailsWithInvalidMutationRate) {
    config.setMutationRate(1.5);
    EXPECT_FALSE(config.validate());
}

//
// 7. Test that getRepresentation returns a non-null pointer.
//
TEST_F(GeneticAlgorithmConfigurationTest, GetRepresentationReturnsNonNull) {
    auto rep = config.getRepresentation();
    EXPECT_NE(rep, nullptr);
}

//
// 8. Test that getSelection returns a non-null pointer.
//
TEST_F(GeneticAlgorithmConfigurationTest, GetSelectionReturnsNonNull) {
    auto sel = config.getSelection();
    EXPECT_NE(sel, nullptr);
}

//
// 9. Test that getCrossover returns a non-null pointer when given an empty generation string.
//
TEST_F(GeneticAlgorithmConfigurationTest, GetCrossoverReturnsNonNull) {
    auto cross = config.getCrossover("");
    EXPECT_NE(cross, nullptr);
}

//
// 10. Test that getMutation returns a non-null pointer.
//
TEST_F(GeneticAlgorithmConfigurationTest, GetMutationReturnsNonNull) {
    auto mut = config.getMutation();
    EXPECT_NE(mut, nullptr);
}

//
// 11. Test that getReplacement returns a non-null pointer. (We set a dummy cache first.)
//
TEST_F(GeneticAlgorithmConfigurationTest, GetReplacementReturnsNonNull) {
    config.setCache(std::make_shared<PopulationCache>());
    auto rep = config.getReplacement();
    EXPECT_NE(rep, nullptr);
}

//
// 12. Test that getStopping returns a non-null pointer.
//
TEST_F(GeneticAlgorithmConfigurationTest, GetStoppingReturnsNonNull) {
    auto stop = config.getStopping();
    EXPECT_NE(stop, nullptr);
}

//
// 13. Test setParameters from a ParameterSet.
// (Assume ParameterSet has an API where you can set values and then later get them as int/double.)
TEST_F(GeneticAlgorithmConfigurationTest, SetParametersFromParameterSet) {
    ParameterSet ps;
    ps.params["populationSize"] = "250";
    ps.params["mutationRate"] = "0.25";
    ps.params["crossoverRate"] = "0.85";
    ps.params["tournamentSize"] = "4";
    config.setParameters(ps);
    EXPECT_EQ(config.getPopulationSize(), 250);
    EXPECT_NEAR(config.getMutationRate(), 0.25, 1e-6);
    EXPECT_NEAR(config.getCrossoverProbability(), 0.85, 1e-6);
    EXPECT_EQ(config.getTournamentSize(), 4);
}

//
// 14. Test setting instance-specific parameters.
//
TEST_F(GeneticAlgorithmConfigurationTest, SetInstanceSpecificParameters) {
    config.setK(10);
    config.setDeltaK(2);
    config.setLNeg(5);
    config.setLPoz(5);
    config.setRepAllowed(false);
    config.setProbablePositive(3);
    EXPECT_EQ(config.getK(), 10);
    EXPECT_EQ(config.getDeltaK(), 2);
    EXPECT_EQ(config.getLNeg(), 5);
    EXPECT_EQ(config.getLPoz(), 5);
    EXPECT_FALSE(config.isRepAllowed());
    EXPECT_EQ(config.getProbablePositive(), 3);
}

//
// 15. Test that adaptive parameters are loaded correctly from a config file.
//
TEST_F(GeneticAlgorithmConfigurationTest, AdaptiveParametersLoadCorrectly) {
    std::string configText =
        "useAdaptiveMutation = false\n"
        "minMutationRate = 0.2\n"
        "maxMutationRate = 0.5\n"
        "stagnationGenerations = 8\n"
        "improvementThreshold = 0.03\n";
    std::string filename = createTempConfigFile(configText);
    bool loaded = config.loadFromFile(filename);
    EXPECT_TRUE(loaded);
    EXPECT_FALSE(config.getAdaptiveParams().useAdaptiveMutation);
    EXPECT_NEAR(config.getAdaptiveParams().minMutationRate, 0.2, 1e-6);
    EXPECT_NEAR(config.getAdaptiveParams().maxMutationRate, 0.5, 1e-6);
    EXPECT_EQ(config.getAdaptiveParams().stagnationGenerations, 8);
    EXPECT_NEAR(config.getAdaptiveParams().improvementThreshold, 0.03, 1e-6);
    removeTempFile(filename);
}

//
// 16. Test that diversity parameters are loaded correctly from a config file.
//
TEST_F(GeneticAlgorithmConfigurationTest, DiversityParametersLoadCorrectly) {
    std::string configText =
        "useFitnessSharing = false\n"
        "useCrowding = true\n"
        "sharingRadius = 0.3\n"
        "sharingAlpha = 1.5\n"
        "diversityWeight = 0.7\n";
    std::string filename = createTempConfigFile(configText);
    bool loaded = config.loadFromFile(filename);
    EXPECT_TRUE(loaded);
    EXPECT_FALSE(config.getDiversityParams().useFitnessSharing);
    EXPECT_TRUE(config.getDiversityParams().useCrowding);
    EXPECT_NEAR(config.getDiversityParams().sharingRadius, 0.3, 1e-6);
    EXPECT_NEAR(config.getDiversityParams().sharingAlpha, 1.5, 1e-6);
    EXPECT_NEAR(config.getDiversityParams().diversityWeight, 0.7, 1e-6);
    removeTempFile(filename);
}

//
// 17. Test that setting the crossover type to "adaptive" causes getCrossover to return an adaptive crossover.
//
TEST_F(GeneticAlgorithmConfigurationTest, CrossoverTypeSettingAdaptive) {
    std::string configText = "crossoverType = adaptive\n";
    std::string filename = createTempConfigFile(configText);
    bool loaded = config.loadFromFile(filename);
    EXPECT_TRUE(loaded);
    auto cross = config.getCrossover("");
    // Here we assume that a dynamic_cast to AdaptiveCrossover is possible.
    auto adaptiveCross = std::dynamic_pointer_cast<AdaptiveCrossover>(cross);
    EXPECT_NE(adaptiveCross, nullptr);
    removeTempFile(filename);
}

//
// 18. Test that an empty selection method causes validate() to fail.
//
TEST_F(GeneticAlgorithmConfigurationTest, ValidateFailsWithEmptySelectionMethod) {
    // Directly set the selection method to empty string.
    config.setSelectionMethod("");
    EXPECT_FALSE(config.validate());
}

//
// 19. Test that the loadFromFile function trims extra whitespace from keys and values.
//
TEST_F(GeneticAlgorithmConfigurationTest, LoadConfigTrimsWhitespace) {
    std::string configText =
        "  populationSize   =    300   \n"
        "mutationRate= 0.3  \n";
    std::string filename = createTempConfigFile(configText);
    bool loaded = config.loadFromFile(filename);
    EXPECT_TRUE(loaded);
    EXPECT_EQ(config.getPopulationSize(), 300);
    EXPECT_NEAR(config.getMutationRate(), 0.3, 1e-6);
    removeTempFile(filename);
}

//
// 20. Test that a malformed line in the config file (without a comma or equals sign) is ignored without crashing.
//
TEST_F(GeneticAlgorithmConfigurationTest, LoadConfigHandlesMalformedLinesGracefully) {
    std::string configText =
        "populationSize=400\n"
        "this is a bad line\n"
        "mutationRate=0.35\n";
    std::string filename = createTempConfigFile(configText);
    bool loaded = config.loadFromFile(filename);
    // Even if there is a malformed line, we expect loadFromFile to complete.
    EXPECT_TRUE(loaded);
    EXPECT_EQ(config.getPopulationSize(), 400);
    EXPECT_NEAR(config.getMutationRate(), 0.35, 1e-6);
    removeTempFile(filename);
}
