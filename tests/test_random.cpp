TEST_CASE("Random number generator with size_t", "[random]") {
    auto& rng = Random::instance();

    SECTION("Basic range tests") {
        // Test small ranges
        for (int i = 0; i < 1000; ++i) {
            size_t value = rng.getRandomSizeT(0, 10);
            REQUIRE(value >= 0);
            REQUIRE(value <= 10);
        }

        // Test larger ranges
        for (int i = 0; i < 1000; ++i) {
            size_t value = rng.getRandomSizeT(100, 1000);
            REQUIRE(value >= 100);
            REQUIRE(value <= 1000);
        }
    }

    SECTION("Edge cases") {
        // Test same min and max
        for (int i = 0; i < 100; ++i) {
            size_t value = rng.getRandomSizeT(5, 5);
            REQUIRE(value == 5);
        }

        // Test consecutive numbers
        for (int i = 0; i < 1000; ++i) {
            size_t value = rng.getRandomSizeT(42, 43);
            REQUIRE(value >= 42);
            REQUIRE(value <= 43);
        }
    }

    SECTION("Distribution tests") {
        // Test distribution for range 0-9
        std::vector<size_t> counts(10, 0);
        const size_t iterations = 10000;
        
        for (size_t i = 0; i < iterations; ++i) {
            size_t value = rng.getRandomSizeT(0, 9);
            REQUIRE(value < counts.size());
            counts[value]++;
        }

        // Check that each number appears at least 8% of the time
        // (allowing for some random variation)
        const size_t minExpected = iterations / 12;  // ~8.3%
        for (size_t i = 0; i < counts.size(); ++i) {
            REQUIRE(counts[i] >= minExpected);
        }
    }

    SECTION("Large number tests") {
        // Test with large numbers that would cause issues if using int
        const size_t largeMin = std::numeric_limits<int>::max();
        const size_t largeMax = largeMin + 1000;

        for (int i = 0; i < 1000; ++i) {
            size_t value = rng.getRandomSizeT(largeMin, largeMax);
            REQUIRE(value >= largeMin);
            REQUIRE(value <= largeMax);
        }
    }

    SECTION("Error cases") {
        // Test max < min
        REQUIRE_THROWS_AS(rng.getRandomSizeT(10, 5), std::invalid_argument);

        // Test with very large ranges
        const size_t maxValue = std::numeric_limits<size_t>::max();
        REQUIRE_NOTHROW(rng.getRandomSizeT(maxValue - 10, maxValue));
    }

    SECTION("Sequence uniqueness test") {
        // Test that we get different numbers in sequence
        const size_t range = 1000;
        std::set<size_t> uniqueValues;
        
        for (int i = 0; i < 100; ++i) {
            size_t value = rng.getRandomSizeT(0, range);
            uniqueValues.insert(value);
        }

        // We should get at least 50 different values in 100 tries
        // (this is a probabilistic test, but the chance of failure is extremely low)
        REQUIRE(uniqueValues.size() > 50);
    }
} 