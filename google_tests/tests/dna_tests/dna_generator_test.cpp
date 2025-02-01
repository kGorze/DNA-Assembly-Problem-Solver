//
// Created by konrad_guest on 01/02/2025.
//
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <future>
#include <filesystem>
#include <sstream>
#include <algorithm>
#include <iostream>
#include "../base_test.h"

// Uwaga: Nie dołączamy dwóch plików definiujących klasę SpectrumGenerator – wystarczy "generator/dna_generator.h".
#include "generator/dna_generator.h"
#include "dna/dna_instance.h"
#include "dna/dna_instance_io.h"
#include "utils/random.h"
#include "utils/logging.h"

// -----------------------------------------------------------------------------------
// Pomocnicza klasa ConstantRandom – zwraca zawsze tą samą wartość (bez override, aby uniknąć błędów)
// -----------------------------------------------------------------------------------
class DummyErrorStrategy : public IErrorIntroductionStrategy {
public:
    void introduceErrors(DNAInstance& instance) override {
        auto spec = instance.getSpectrum();
        if (!spec.empty()) {
            spec.erase(spec.begin());
            instance.setSpectrum(spec);
        }
    }
};

class ConstantRandom : public Random {
public:
    ConstantRandom() : Random() {}
    explicit ConstantRandom(int seed) : Random(seed) {}
    
    // Always returns 1 (nucleotides[1]=='C')
    int getRandomInt(int /*min*/, int /*max*/) override { return 1; }
    
    // Always returns 0.5 for probability
    double generateProbability() override { return 0.5; }
    
    // Always returns 1 for size_t as well
    size_t getRandomSizeT(size_t /*min*/, size_t /*max*/) override { return 1; }
};

// -----------------------------------------------------------------------------------
// Test Fixture for DNAGenerator tests - handles logger initialization/cleanup
// -----------------------------------------------------------------------------------
class DNAGeneratorTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("dna_generator_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();  // Call base class SetUp
    }

    void TearDown() override {
        BaseTest::TearDown();  // Call base class TearDown
        Logger::cleanup();
    }
};

// -----------------------------------------------------------------------------------
// Testy DNAGenerator (25 testów)
// -----------------------------------------------------------------------------------

// 1. Konstruktor – tworzymy DNAGenerator z unikalnym Random.
TEST_F(DNAGeneratorTest, ConstructorWithValidRandom) {
    EXPECT_NO_THROW({
        DNAGenerator gen(std::make_unique<Random>());
    });
}

// 2. Konstruktor z nullptr – powinien utworzyć domyślny generator
TEST_F(DNAGeneratorTest, ConstructorWithNullRandom) {
    DNAGenerator gen(nullptr);
    std::string dna = gen.generateDNA(10, true);
    EXPECT_EQ(dna.length(), 10);
}

// 3. setParameters z poprawnymi wartościami – nie rzuca wyjątku.
TEST_F(DNAGeneratorTest, SetParametersValid) {
    DNAGenerator gen(std::make_unique<Random>());
    
    EXPECT_NO_THROW({
         gen.setParameters(20, 5, 1);
    });
    
    // Aby sprawdzić poprawność, wywołamy generateRandomInstance – powinno działać.
    EXPECT_NO_THROW({
         DNAInstance inst = gen.generateRandomInstance(20, 5, 2, 1);
         EXPECT_EQ(inst.getDNA().length(), 20);
         
         auto spectrum = inst.getSpectrum();
         for (const auto& kmer : spectrum) {
             LOG_INFO("k-mer: {}", kmer);
         }
    });
}

// 4. setParameters: n <= 0 – rzuca wyjątek.
TEST_F(DNAGeneratorTest, SetParametersInvalid_N) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_THROW({
         gen.setParameters(0, 10, 2);
    }, std::invalid_argument);
}

// 5. setParameters: k <= 0 – rzuca wyjątek.
TEST_F(DNAGeneratorTest, SetParametersInvalid_K) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_THROW({
         gen.setParameters(100, 0, 2);
    }, std::invalid_argument);
}

// 6. setParameters: deltaK < 0 – rzuca wyjątek.
TEST_F(DNAGeneratorTest, SetParametersInvalid_DeltaK) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_THROW({
         gen.setParameters(100, 10, -1);
    }, std::invalid_argument);
}

// 7. setParameters: k > n – rzuca wyjątek.
TEST_F(DNAGeneratorTest, SetParametersInvalid_KGreaterThanN) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_THROW({
         gen.setParameters(50, 60, 2);
    }, std::invalid_argument);
}

// 8. generateDNA z poprawną długością – wynikowy ciąg ma zadaną długość.
TEST_F(DNAGeneratorTest, GenerateDNAValidLength) {
    DNAGenerator gen(std::make_unique<Random>());
    std::string dna = gen.generateDNA(20, true);
    EXPECT_EQ(dna.length(), 20);
}

// 9. generateDNA: długość <= 0 – rzuca wyjątek.
TEST_F(DNAGeneratorTest, GenerateDNAInvalidLength) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_THROW({
         gen.generateDNA(0, true);
    }, std::invalid_argument);
}

// 10. generateDNA: wygenerowany ciąg zawiera tylko dozwolone znaki.
TEST_F(DNAGeneratorTest, GenerateDNAValidCharacters) {
    DNAGenerator gen(std::make_unique<Random>());
    std::string dna = gen.generateDNA(50, true);
    const std::string validChars = "ACGT";
    for (char c : dna) {
         EXPECT_NE(validChars.find(c), std::string::npos)
             << "Nieprawidłowy znak: " << c;
    }
}

// 11. generateDNA z repAllowed == false – mimo braku implementacji repAllowed, ciąg jest poprawny.
TEST_F(DNAGeneratorTest, GenerateDNARepetitionNotAllowed) {
    DNAGenerator gen(std::make_unique<Random>());
    std::string dna = gen.generateDNA(30, false);
    EXPECT_EQ(dna.length(), 30);
    const std::string validChars = "ACGT";
    for (char c : dna) {
         EXPECT_NE(validChars.find(c), std::string::npos);
    }
}

// 12. Indirect test: po ustawieniu prawidłowych parametrów, metoda validateParameters (prywatna) działa – poprzez generateRandomInstance.
TEST_F(DNAGeneratorTest, GenerateRandomInstanceValidParameters) {
    DNAGenerator gen(std::make_unique<Random>());
    gen.setParameters(120, 12, 2);
    EXPECT_NO_THROW({
         DNAInstance inst = gen.generateRandomInstance(120, 12, 4, 4);
         EXPECT_EQ(inst.getDNA().length(), 120);
    });
}

// 13. generateRandomInstance: wynikowy DNAInstance ma długość równą size.
TEST_F(DNAGeneratorTest, GenerateRandomInstanceCorrectSize) {
    DNAGenerator gen(std::make_unique<Random>());
    int size = 150, k = 10;
    gen.setParameters(size, k, 3);
    DNAInstance inst = gen.generateRandomInstance(size, k, 5, 5);
    EXPECT_EQ(inst.getDNA().length(), size);
}

// 14. generateRandomInstance: spectrum powinno mieć (size - k + 1) k-merów.
TEST_F(DNAGeneratorTest, GenerateRandomInstanceSpectrumSize) {
    DNAGenerator gen(std::make_unique<Random>());
    int size = 150, k = 10;
    gen.setParameters(size, k, 3);
    DNAInstance inst = gen.generateRandomInstance(size, k, 5, 5);
    EXPECT_EQ(inst.getSpectrum().size(), size - k + 1);
}

// 15. saveToFile – przy poprawnych danych metoda zwraca true.
TEST_F(DNAGeneratorTest, SaveToFileSuccess) {
    DNAGenerator gen(std::make_unique<Random>());
    DNAInstance instance(100, 10, 3, 2, 0, false, 0.0, 0);
    bool res = gen.saveToFile(instance, "temp_instance.txt");
    EXPECT_TRUE(res);
    std::filesystem::remove("temp_instance.txt");
}

// 16. saveToFile – przy pustym filename metoda zwraca false (lub nie zapisuje).
TEST_F(DNAGeneratorTest, SaveToFileFailureEmptyFilename) {
    DNAGenerator gen(std::make_unique<Random>());
    DNAInstance instance(100, 10, 3, 2, 0, false, 0.0, 0);
    bool res = gen.saveToFile(instance, "");
    EXPECT_FALSE(res);
}

// 17. loadFromFile – zapisujemy instancję do pliku, potem ją odczytujemy.
TEST_F(DNAGeneratorTest, LoadFromFileSuccess) {
    DNAGenerator gen(std::make_unique<Random>());
    DNAInstance instance(100, 10, 3, 2, 0, false, 0.0, 0);
    instance.setDNA("ACGTACGTAC");
    std::vector<std::string> spec = {"ACGTACGTAC", "CGTACGTACG", "GTACGTACGT", "TACGTACGTA", "ACGTACGTAC"};
    instance.setSpectrum(spec);
    bool saved = gen.saveToFile(instance, "temp_instance.txt");
    EXPECT_TRUE(saved);
    DNAInstance loaded = DNAGenerator::loadFromFile("temp_instance.txt");
    EXPECT_EQ(loaded.getDNA(), instance.getDNA());
    std::filesystem::remove("temp_instance.txt");
}

// 18. loadFromFile – dla nieistniejącego pliku metoda rzuca wyjątek.
TEST_F(DNAGeneratorTest, LoadFromFileFailureNonexistentFile) {
    EXPECT_THROW({
         DNAGenerator::loadFromFile("nonexistent_file.txt");
    }, std::runtime_error);
}

// 19. generateDNASpectrum – dla zadanej instancji metoda zwraca oczekiwane k-mery.
TEST_F(DNAGeneratorTest, GenerateDNASpectrumValid) {
    // Ustalony DNA: "ACGTAC"
    // k = 3, deltaK = 0
    DNAInstance instance(6, 3, 0, 0, 0, false, 0.0, 0);
    instance.setDNA("ACGTAC");
    // Oczekiwane k-mery: "ACG", "CGT", "GTA", "TAC"
    DNAGenerator gen(std::make_unique<Random>());
    auto spectrum = gen.generateDNASpectrum(instance);
    std::vector<std::string> expected = {"ACG", "CGT", "GTA", "TAC"};
    EXPECT_EQ(spectrum, expected);
}

// 20. generateDNASpectrum – jeżeli DNA jest puste, metoda rzuci wyjątek.
TEST_F(DNAGeneratorTest, GenerateDNASpectrumEmptyDNA) {
    // Tworzymy nową instancję - domyślnie DNA będzie puste
    DNAInstance instance;
    DNAGenerator gen(std::make_unique<Random>());
    
    // Sprawdzamy czy generateDNASpectrum rzuci wyjątek dla instancji z pustym DNA
    EXPECT_THROW({
        gen.generateDNASpectrum(instance);
    }, std::invalid_argument);
}

// 21. Test współbieżności generateDNA – wywołanie metody z kilku wątków.
TEST_F(DNAGeneratorTest, GenerateDNAThreadSafety) {
    DNAGenerator gen(std::make_unique<Random>());
    const int length = 50;
    const int numThreads = 4;
    std::vector<std::future<std::string>> futures;
    for (int i = 0; i < numThreads; ++i) {
        futures.push_back(std::async(std::launch::async, [&gen, length]() {
            return gen.generateDNA(length, true);
        }));
    }
    for (auto &fut : futures) {
        std::string dna = fut.get();
        EXPECT_EQ(dna.length(), length);
    }
}

// 22. Test deterministyczności – korzystając z ConstantRandom, wygenerowany DNA powinien być taki sam.
TEST_F(DNAGeneratorTest, GenerateDNADeterministic) {
    auto constantRand = std::make_unique<ConstantRandom>();
    DNAGenerator gen(std::move(constantRand));
    std::string dna = gen.generateDNA(20, true);
    std::string expected(20, 'C'); // nucleotides[1] to 'C'
    EXPECT_EQ(dna, expected);
}

// 23. Test – jeżeli parametry nie zostały ustawione (validateParameters false), generateRandomInstance rzuca wyjątek.
TEST_F(DNAGeneratorTest, GenerateRandomInstanceWithoutSettingParameters) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_THROW({
         gen.generateRandomInstance(150, 10, 5, 5);
    }, std::invalid_argument);
}

// 24. Test generateRandomInstance – sprawdzenie, czy spectrum instancji ma poprawne długości k-merów.
TEST_F(DNAGeneratorTest, GenerateRandomInstanceSpectrumLength) {
    DNAGenerator gen(std::make_unique<Random>());
    int size = 150, k = 10;
    gen.setParameters(size, k, 3);
    DNAInstance inst = gen.generateRandomInstance(size, k, 5, 5);
    auto spectrum = inst.getSpectrum();
    for (const auto &oligo : spectrum) {
         EXPECT_EQ(oligo.length(), k);
    }
}

// 25. Dodatkowy test: sprawdzenie, że generateDNA nie powoduje awarii przy repAllowed == false.
TEST_F(DNAGeneratorTest, GenerateDNARepNotAllowedDoesNotCrash) {
    DNAGenerator gen(std::make_unique<Random>());
    EXPECT_NO_THROW({
         std::string dna = gen.generateDNA(40, false);
         EXPECT_EQ(dna.length(), 40);
    });
}
