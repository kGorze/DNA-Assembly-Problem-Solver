//
// Created by konrad_guest on 01/02/2025.
//
// tests/utils/random_test.cpp
#include <gtest/gtest.h>
#include "utils/random.h"  // Upewnij się, że ścieżka do pliku nagłówkowego jest poprawna
#include "../base_test.h"
#include <future>
#include <thread>
#include <vector>
#include <stdexcept>

// Przyjmujemy, że klasa Random posiada konstruktor domyślny, konstruktor przyjmujący seed,
// metody: getRandomInt(int min, int max), generateProbability(), setSeed(unsigned int)
// oraz (opcjonalnie) metodę getRandomSizeT(min, max).

class RandomTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("random_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();
    }

    void TearDown() override {
        BaseTest::TearDown();
        Logger::cleanup();
    }
};

// Test 1: Sprawdzenie, że getRandomInt zwraca wartość z zadanego przedziału.
TEST_F(RandomTest, GetRandomIntInRange) {
    Random r; // konstruktor domyślny
    for (int i = 0; i < 1000; i++) {
        int val = r.getRandomInt(5, 15);
        EXPECT_GE(val, 5);
        EXPECT_LE(val, 15);
    }
}

// Test 2: Sprawdzenie, że generateProbability zwraca wartość z przedziału [0,1].
TEST_F(RandomTest, GenerateProbabilityRange) {
    Random r;
    for (int i = 0; i < 1000; i++) {
        double p = r.generateProbability();
        EXPECT_GE(p, 0.0);
        EXPECT_LE(p, 1.0);
    }
}

// Test 3: Powtarzalność – dwa obiekty Random z tym samym seed'em generują tę samą sekwencję.
TEST_F(RandomTest, Reproducibility) {
    unsigned int seed = 12345;
    Random r1(seed);
    Random r2(seed);
    const int iterations = 100;
    for (int i = 0; i < iterations; i++) {
        int val1 = r1.getRandomInt(0, 100);
        int val2 = r2.getRandomInt(0, 100);
        EXPECT_EQ(val1, val2);
    }
}

// Test 4: Test that two Random instances with the same seed generate identical sequences
TEST_F(RandomTest, SameSeedGeneratesSameSequence) {
    // Create two Random instances with the same seed
    unsigned int seed = 98765;
    Random r1(seed);
    Random r2(seed);
    
    // Generate sequences from both instances
    std::vector<int> seq1, seq2;
    for (int i = 0; i < 50; i++) {
        seq1.push_back(r1.getRandomInt(0, 50));
        seq2.push_back(r2.getRandomInt(0, 50));
    }
    
    // The sequences should be identical
    EXPECT_EQ(seq1, seq2);
}

// Test 5: Test współbieżności – wywołanie getRandomInt w kilku wątkach.
TEST_F(RandomTest, ThreadSafety) {
    Random r;
    const int iterations = 100;
    const int numThreads = 4;
    std::vector<std::future<std::vector<int>>> futures;
    for (int t = 0; t < numThreads; ++t) {
        futures.push_back(std::async(std::launch::async, [&r, iterations]() {
            std::vector<int> results;
            for (int i = 0; i < iterations; i++) {
                results.push_back(r.getRandomInt(10, 20));
            }
            return results;
        }));
    }
    for (auto &f : futures) {
        auto vec = f.get();
        for (int val : vec) {
            EXPECT_GE(val, 10);
            EXPECT_LE(val, 20);
        }
    }
}

// Test 6: Jeśli min == max, getRandomInt powinien zwrócić tę wartość.
TEST_F(RandomTest, GetRandomIntSameMinMax) {
    Random r;
    for (int i = 0; i < 100; i++) {
        int val = r.getRandomInt(7, 7);
        EXPECT_EQ(val, 7);
    }
}

// Test 7: Jeśli metoda getRandomInt wywołana z nieprawidłowym przedziałem (min > max) powinna rzucać wyjątek.
// (Jeśli Twoja implementacja tego nie robi, możesz pominąć lub zmodyfikować ten test.)
TEST_F(RandomTest, GetRandomIntInvalidRange) {
    Random r;
    EXPECT_THROW({
        r.getRandomInt(10, 5);
    }, std::invalid_argument);
}

// Test 8: Sprawdzamy, że metoda generateProbability nie zwraca stale tej samej wartości.
TEST_F(RandomTest, GenerateProbabilityVariability) {
    Random r;
    std::vector<double> probs;
    for (int i = 0; i < 1000; i++) {
        probs.push_back(r.generateProbability());
    }
    double first = probs.front();
    bool allSame = std::all_of(probs.begin(), probs.end(), [first](double v) { return v == first; });
    EXPECT_FALSE(allSame);
}

// Test 9: Test, że domyślny konstruktor (bez zadanego seed'u) nie generuje deterministycznej sekwencji.
TEST_F(RandomTest, DefaultConstructorNonDeterminism) {
    Random r1;
    Random r2;
    bool allEqual = true;
    for (int i = 0; i < 100; i++) {
        if (r1.getRandomInt(0, 100) != r2.getRandomInt(0, 100)) {
            allEqual = false;
            break;
        }
    }
    EXPECT_FALSE(allEqual);
}

// Test 10: Jeśli klasa Random posiada metodę getRandomSizeT (podobną do getRandomInt dla size_t),
// sprawdzamy, czy zwraca wartości w zadanym przedziale.
TEST_F(RandomTest, GetRandomSizeTInRange) {
    Random r;
    for (int i = 0; i < 1000; i++) {
        size_t val = r.getRandomSizeT(10, 20);
        EXPECT_GE(val, 10);
        EXPECT_LE(val, 20);
    }
}
