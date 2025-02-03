//
// Created by konrad_guest on 02/02/2025. - test already in the suite
//

#include <gtest/gtest.h>
#include "../base_test.h"           // Używamy BaseTest jako fixture
#include "../../../include/configuration/genetic_algorithm_configuration.h"
#include "../../../include/metaheuristics/adaptive_mutation.h"        // Plik z definicją AdaptiveMutation
#include "../../../include/utils/logging.h"  // Poprawiona ścieżka do logging.h
#include <memory>

// Forward declaration
struct AdaptiveParams;

// DummyGAConfig – uproszczona konfiguracja, wystarczająca do testowania AdaptiveMutation
class DummyGAConfig final : public GAConfig {
public:
    // Factory method to ensure proper shared_ptr creation
    static std::shared_ptr<DummyGAConfig> create(double mutationRate, const AdaptiveParams& params) {
        return std::shared_ptr<DummyGAConfig>(new DummyGAConfig(mutationRate, params));
    }

private:
    // Make constructor private to force use of factory method
    explicit DummyGAConfig(double mutationRate, const AdaptiveParams& params) : GAConfig() {
        LOG_DEBUG("DummyGAConfig constructor - creating with mutation rate: {}", mutationRate);
        setMutationRate(mutationRate);
        setAdaptiveParams(params);
    }

    // Removed friend declarations since factory method now uses new directly
    // friend class std::allocator<DummyGAConfig>;
    // template<typename _Alloc>
    // friend class std::_Sp_counted_ptr_inplace;

public:
    ~DummyGAConfig() override {
        LOG_DEBUG("DummyGAConfig destructor - destroying at address: {}", static_cast<void*>(this));
    }
    
    // Delete copy constructor and assignment to prevent accidental copies
    DummyGAConfig(const DummyGAConfig&) = delete;
    DummyGAConfig& operator=(const DummyGAConfig&) = delete;
    DummyGAConfig(DummyGAConfig&&) = delete;
    DummyGAConfig& operator=(DummyGAConfig&&) = delete;
};

// Test fixture – wykorzystujemy BaseTest, aby uzyskać standardowy output z testów.
class AdaptiveMutationTest : public BaseTest {
protected:
    static bool s_loggerInitialized;
    
    // Move member variables to protected section
    std::shared_ptr<DummyGAConfig> config;
    std::unique_ptr<AdaptiveMutation> adaptive;

    static void SetUpTestSuite() {
        if (!s_loggerInitialized) {
            try {
                Logger::initialize("adaptive_mutation_test.log");
                Logger::setLogLevel(LogLevel::DEBUG);
                s_loggerInitialized = true;
                LOG_INFO("AdaptiveMutationTest test suite initialized successfully");
            } catch (const std::exception& e) {
                std::cerr << "Failed to initialize logger: " << e.what() << std::endl;
                throw;
            }
        }
    }

    static void TearDownTestSuite() {
        if (s_loggerInitialized) {
            LOG_INFO("Tearing down AdaptiveMutationTest test suite");
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
            // Initialize member variables to nullptr first
            config = nullptr;
            adaptive = nullptr;
            
            // Ensure logger is initialized
            if (!s_loggerInitialized) {
                SetUpTestSuite();
            }
            
            LOG_DEBUG("AdaptiveMutationTest::SetUp - Starting setup");
            BaseTest::SetUp();
            
            AdaptiveParams params{};  
            params.useAdaptiveMutation = true;
            params.improvementThreshold = 1.0;
            params.minMutationRate = 0.01;
            params.maxMutationRate = 0.5;
            params.stagnationGenerations = 3;
            
            LOG_DEBUG("Creating DummyGAConfig");
            // Use factory method to create config
            config = DummyGAConfig::create(0.1, params);
            if (!config) {
                LOG_ERROR("Failed to create DummyGAConfig - nullptr returned");
                throw std::runtime_error("Failed to create DummyGAConfig");
            }
            LOG_DEBUG("DummyGAConfig created successfully at address: {}", static_cast<void*>(config.get()));
            
            LOG_DEBUG("Creating AdaptiveMutation");
            adaptive = std::make_unique<AdaptiveMutation>(config);
            if (!adaptive) {
                LOG_ERROR("Failed to create AdaptiveMutation - nullptr returned");
                throw std::runtime_error("Failed to create AdaptiveMutation");
            }
            LOG_DEBUG("AdaptiveMutation created successfully at address: {}", static_cast<void*>(adaptive.get()));
            
            LOG_DEBUG("AdaptiveMutationTest::SetUp - Setup completed successfully");
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
            LOG_DEBUG("AdaptiveMutationTest::TearDown - Starting cleanup");
            
            // (1) First clear any references in AdaptiveMutation
            if (adaptive) {
                LOG_DEBUG("Clearing references in AdaptiveMutation at address: {}", static_cast<void*>(adaptive.get()));
                adaptive->clearReferences();
            }
            
            // (2) Then destroy AdaptiveMutation
            if (adaptive) {
                LOG_DEBUG("Releasing AdaptiveMutation at address: {}", static_cast<void*>(adaptive.get()));
                adaptive.reset();  // Use reset directly instead of move+reset
            }
            
            // (3) Now safe to destroy config since AdaptiveMutation is gone
            if (config) {
                LOG_DEBUG("Releasing DummyGAConfig at address: {}", static_cast<void*>(config.get()));
                config.reset();  // Use reset directly instead of move+reset
            }
            
            LOG_DEBUG("Calling BaseTest::TearDown");
            BaseTest::TearDown();
            LOG_DEBUG("AdaptiveMutationTest::TearDown - Cleanup completed successfully");
        } catch (const std::exception& e) {
            LOG_ERROR("TearDown failed with exception: {}", e.what());
            std::cerr << "TearDown failed with exception: " << e.what() << std::endl;
            // Don't rethrow, just log the error
            ADD_FAILURE() << "TearDown failed: " << e.what();
        } catch (...) {
            LOG_ERROR("TearDown failed with unknown exception");
            std::cerr << "TearDown failed with unknown exception" << std::endl;
            // Don't rethrow, just log the error
            ADD_FAILURE() << "TearDown failed with unknown exception";
        }
    }
    
    void recreateObjects(double mutationRate, const AdaptiveParams& params) {
        try {
            LOG_DEBUG("recreateObjects - Starting recreation with mutation rate: {}", mutationRate);
            
            // (1) First clear any references in AdaptiveMutation
            if (adaptive) {
                LOG_DEBUG("Clearing references in AdaptiveMutation at address: {}", static_cast<void*>(adaptive.get()));
                adaptive->clearReferences();
            }
            
            // (2) Then destroy AdaptiveMutation
            if (adaptive) {
                LOG_DEBUG("Releasing old AdaptiveMutation at address: {}", static_cast<void*>(adaptive.get()));
                adaptive.reset();  // Use reset directly instead of move+reset
            }
            
            // (3) Then destroy config
            if (config) {
                LOG_DEBUG("Releasing old DummyGAConfig at address: {}", static_cast<void*>(config.get()));
                config.reset();  // Use reset directly instead of move+reset
            }
            
            // (4) Create fresh config first using factory method
            LOG_DEBUG("Creating new DummyGAConfig");
            config = DummyGAConfig::create(mutationRate, params);
            if (!config) {
                LOG_ERROR("Failed to recreate DummyGAConfig - nullptr returned");
                throw std::runtime_error("Failed to recreate DummyGAConfig");
            }
            LOG_DEBUG("New DummyGAConfig created at address: {}", static_cast<void*>(config.get()));
            
            // (5) Then create AdaptiveMutation that depends on config
            LOG_DEBUG("Creating new AdaptiveMutation");
            adaptive = std::make_unique<AdaptiveMutation>(config);
            if (!adaptive) {
                LOG_ERROR("Failed to recreate AdaptiveMutation - nullptr returned");
                throw std::runtime_error("Failed to recreate AdaptiveMutation");
            }
            LOG_DEBUG("New AdaptiveMutation created at address: {}", static_cast<void*>(adaptive.get()));
            
            LOG_DEBUG("recreateObjects - Recreation completed successfully");
        } catch (const std::exception& e) {
            LOG_ERROR("recreateObjects failed with exception: {}", e.what());
            FAIL() << "recreateObjects failed: " << e.what();
        } catch (...) {
            LOG_ERROR("recreateObjects failed with unknown exception");
            FAIL() << "recreateObjects failed with unknown exception";
        }
    }

    // Gettery zwracają referencje do wskaźników
    const std::shared_ptr<DummyGAConfig>& getConfig() const { 
        LOG_DEBUG("Getting config at address: {}", static_cast<void*>(config.get()));
        return config; 
    }
    
    const std::unique_ptr<AdaptiveMutation>& getAdaptive() const { 
        LOG_DEBUG("Getting adaptive at address: {}", static_cast<void*>(adaptive.get()));
        return adaptive; 
    }
};

// Initialize static member
bool AdaptiveMutationTest::s_loggerInitialized = false;

//
// Test 1: Początkowa wartość mutacji musi być równa tej z konfiguracji
//
TEST_F(AdaptiveMutationTest, InitialMutationRate) {
    try {
        LOG_DEBUG("Starting InitialMutationRate test");
        
        LOG_DEBUG("Checking if adaptive is not null");
        const auto& adaptive_ref = getAdaptive();
        ASSERT_TRUE(adaptive_ref != nullptr) << "Adaptive mutation object is null";
        LOG_DEBUG("Adaptive mutation object is valid at address: {}", static_cast<void*>(adaptive_ref.get()));
        
        LOG_DEBUG("Checking if config is not null");
        const auto& config_ref = getConfig();
        ASSERT_TRUE(config_ref != nullptr) << "Config object is null";
        LOG_DEBUG("Config object is valid at address: {}", static_cast<void*>(config_ref.get()));
        
        LOG_DEBUG("Getting current mutation rate");
        double current_rate = adaptive_ref->getCurrentMutationRate();
        LOG_DEBUG("Current mutation rate is: {}", current_rate);
        
        LOG_DEBUG("Comparing mutation rates");
        EXPECT_NEAR(current_rate, 0.1, 1e-6);
        
        LOG_DEBUG("InitialMutationRate test completed successfully");
    } catch (const std::exception& e) {
        LOG_ERROR("InitialMutationRate test failed with exception: {}", e.what());
        throw;
    } catch (...) {
        LOG_ERROR("InitialMutationRate test failed with unknown exception");
        throw;
    }
}

//
// Test 2: Przy znaczącej poprawie (improvement > threshold) mutacja maleje – mnożenie przez 0.95
//
TEST_F(AdaptiveMutationTest, ImprovementDecreasesMutationRate) {
    getAdaptive()->reset();
    // Pierwsze wywołanie: bestFitness = 2.0, poprawa = 2.0 - 0.0 = 2.0 > threshold (1.0)
    getAdaptive()->updateMutationRate(2.0);
    double expected = 0.1 * 0.95; // 0.095 (jest większe niż min 0.01)
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), expected, 1e-6);
}

//
// Test 3: Brak poprawy (improvement <= threshold) – po 3 wywołaniach (stagnationGenerations) wartość mutacji wzrasta
//
TEST_F(AdaptiveMutationTest, NoImprovementIncreasesMutationRateAfterStagnation) {
    getAdaptive()->reset();
    // 1. update: bestFitness = 0.5 (poprawa = 0.5, nie > 1.0) → stagnacja = 1, mut rate = 0.1
    getAdaptive()->updateMutationRate(0.5);
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), 0.1, 1e-6);
    // 2. update: bestFitness = 0.6 (poprawa = 0.1) → stagnacja = 2, mut rate = 0.1
    getAdaptive()->updateMutationRate(0.6);
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), 0.1, 1e-6);
    // 3. update: bestFitness = 0.7 (poprawa = 0.1) → stagnacja osiąga 3, mut rate wzrasta: 0.1 * 1.1 = 0.11
    getAdaptive()->updateMutationRate(0.7);
    double expected = 0.1 * 1.1;
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), expected, 1e-6);
}

//
// Test 4: Sprawdzenie, czy updateMutationRate aktualizuje wewnętrznie ostatnią wartość fitness
//
TEST_F(AdaptiveMutationTest, UpdateLastBestFitness) {
    getAdaptive()->reset();
    getAdaptive()->updateMutationRate(3.0);
    // Kolejna aktualizacja z wartością 3.5 (poprawa = 0.5, nie > threshold) → mutation rate nie zmieni się
    getAdaptive()->updateMutationRate(3.5);
    double expected = 0.1 * 0.95;  // Zgodnie z pierwszym update (bo brak zmiany)
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), expected, 1e-6);
}

//
// Test 5: Gdy adaptive mutation jest wyłączone, updateMutationRate nie zmienia wartości mutacji
//
TEST_F(AdaptiveMutationTest, AdaptiveMutationDisabled) {
    // Wyłączamy adaptive mutation poprzez zmianę parametru
    AdaptiveParams params{};
    params.useAdaptiveMutation = false;
    
    // Używamy helper metody do bezpiecznego przetworzenia obiektów
    recreateObjects(0.1, params);
    
    ASSERT_NE(getAdaptive(), nullptr);
    ASSERT_NE(getConfig(), nullptr);
    
    getAdaptive()->updateMutationRate(5.0);
    // Oczekujemy, że mutacja pozostanie niezmieniona
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), 0.1, 1e-6);
}

//
// Test 6: Metoda reset() przywraca stan początkowy
//
TEST_F(AdaptiveMutationTest, ResetResetsState) {
    getAdaptive()->updateMutationRate(2.0);  // Zmiana stanu
    getAdaptive()->reset();
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), 0.1, 1e-6);
    // Po resecie, updateMutationRate z poprawą powinno obniżyć wartość mutacji
    getAdaptive()->updateMutationRate(2.0);
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), 0.1 * 0.95, 1e-6);
}

//
// Test 7: Jeśli poprawa jest równa dokładnie threshold, nie traktujemy tego jako poprawy
//
TEST_F(AdaptiveMutationTest, ImprovementEqualToThresholdNoDecrease) {
    getAdaptive()->reset();
    // Pierwsze update: bestFitness = 1.0, poprawa = 1.0, nie > 1.0
    getAdaptive()->updateMutationRate(1.0);
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), 0.1, 1e-6);
}

//
// Test 8: Wielokrotne sukcesywne poprawy powodują skumulowane zmniejszenie mutation rate
//
TEST_F(AdaptiveMutationTest, MultipleImprovementsCompoundDecrease) {
    getAdaptive()->reset();
    getAdaptive()->updateMutationRate(2.0);  // 0.1 * 0.95 = 0.095
    getAdaptive()->updateMutationRate(4.0);  // kolejne mnożenie: 0.095 * 0.95
    double expected = 0.1 * 0.95 * 0.95;
    EXPECT_NEAR(getAdaptive()->getCurrentMutationRate(), expected, 1e-6);
}

//
// Test 9: Mutation rate nigdy nie spada poniżej minimalnej wartości
//
TEST_F(AdaptiveMutationTest, MutationRateNotBelowMin) {
    getAdaptive()->reset();
    // Symulujemy wiele update'ów z poprawami, które będą zmniejszać rate
    for (int i = 0; i < 50; ++i) {
        getAdaptive()->updateMutationRate(10.0 + i);
    }
    const auto& params = getConfig()->getAdaptiveParams();
    EXPECT_GE(getAdaptive()->getCurrentMutationRate(), params.minMutationRate);
}

//
// Test 10: Mutation rate nigdy nie przekracza maksymalnej wartości przy wzroście po stagnacji
//
TEST_F(AdaptiveMutationTest, MutationRateNotAboveMax) {
    getAdaptive()->reset();
    // Symulujemy wiele update'ów bez poprawy (stagnacja)
    for (int i = 0; i < 10; ++i) {
        getAdaptive()->updateMutationRate(0.5);  // Brak istotnej poprawy
    }
    const auto& params = getConfig()->getAdaptiveParams();
    EXPECT_LE(getAdaptive()->getCurrentMutationRate(), params.maxMutationRate);
}
