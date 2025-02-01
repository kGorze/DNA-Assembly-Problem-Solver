//
// Created by konrad_guest on 01/02/2025.
//
#include <gtest/gtest.h>
#include "dna/spectrum_generator.h"
#include "utils/logging.h"
#include "dna/dna_instance.h"
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>
#include "../base_test.h"  // Zakładamy, że BaseTest ustawia Logger, itp.

class SpectrumGeneratorTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("spectrum_generator_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();
    }
    void TearDown() override {
        BaseTest::TearDown();
        Logger::cleanup();
    }
};

// Test 1: Konstruktor – sprawdzenie, że obiekt tworzy się przy podaniu wskaźnika do deterministycznego generatora.
TEST_F(SpectrumGeneratorTest, ConstructorInitializesRandom) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    EXPECT_NO_THROW({
        SpectrumGenerator sg(rng);
    });
}

// Test 2: generateSpectrum rzuca wyjątek przy pustym DNA.
TEST_F(SpectrumGeneratorTest, GenerateSpectrumEmptyDNAThrows) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    EXPECT_THROW({
        sg.generateSpectrum("", 5, 0);
    }, std::invalid_argument);
}

// Test 3: generateSpectrum rzuca wyjątek, gdy k jest niepoprawne (k <= 0 lub k > długość DNA).
TEST_F(SpectrumGeneratorTest, GenerateSpectrumInvalidKThrows) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    std::string dna = "ACGTACGT";
    EXPECT_THROW({
        sg.generateSpectrum(dna, 0, 0);
    }, std::invalid_argument);
    EXPECT_THROW({
        sg.generateSpectrum(dna, 20, 0);
    }, std::invalid_argument);
}

// Test 4: generateSpectrum rzuca wyjątek, gdy errorCount < 0.
TEST_F(SpectrumGeneratorTest, GenerateSpectrumNegativeErrorCountThrows) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    std::string dna = "ACGTACGT";
    EXPECT_THROW({
        sg.generateSpectrum(dna, 3, -1);
    }, std::invalid_argument);
}

// Test 5: Generowanie spekturm-u bez błędów – sprawdzamy, że metoda zwraca unikalne k-mery.
TEST_F(SpectrumGeneratorTest, GenerateSpectrumWithoutErrors) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    std::string dna = "ACGTACGT"; // długość 8
    int k = 3;
    int errorCount = 0;
    // Oryginalne k-mery: "ACG", "CGT", "GTA", "TAC", "ACG", "CGT"
    // Unikalne: {"ACG", "CGT", "GTA", "TAC"}
    std::vector<std::string> expected = {"ACG", "CGT", "GTA", "TAC"};
    std::sort(expected.begin(), expected.end());
    std::vector<std::string> spectrum = sg.generateSpectrum(dna, k, errorCount);
    std::sort(spectrum.begin(), spectrum.end());
    EXPECT_EQ(spectrum, expected);
}

// Test 6: Generowanie spekturm-u z błędami – sprawdzamy, że do oryginalnych k-merów dodane zostają error warianty.
// Dla krótkiego DNA łatwiej oszacować oczekiwany zbiór.
TEST_F(SpectrumGeneratorTest, GenerateSpectrumWithErrors) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    // Użyjemy DNA o długości 3, k = 2
    std::string dna = "ACG"; // możliwe k-mery: "AC" i "CG"
    int k = 2;
    int errorCount = 1;  // efektywnie generujemy warianty z jedną zamianą
    // Dla "AC": zmiana pozycji 0 → ("CC", "GC", "TC"), pozycji 1 → ("AA", "AG", "AT")
    // Dla "CG": zmiana pozycji 0 → ("AG", "GG", "TG"), pozycji 1 → ("CA", "CC", "CT")
    // Dodajemy oryginalne k-mery "AC" i "CG"
    std::vector<std::string> expected;
    // Zbiór oczekiwany może zawierać więcej pozycji – zamiast dokładnie sprawdzać każdy wariant,
    // sprawdzimy, że wynikowy spekturm zawiera oryginały oraz ma więcej elementów niż bez błędów.
    std::vector<std::string> spectrumWithoutErrors = sg.generateSpectrum(dna, k, 0);
    std::vector<std::string> spectrumWithErrors = sg.generateSpectrum(dna, k, errorCount);
    // Spektrum z błędami powinno być rozszerzone względem spektrum bez błędów
    EXPECT_GT(spectrumWithErrors.size(), spectrumWithoutErrors.size());
    // Sprawdzamy, czy oryginalne k-mery znajdują się w wyniku
    for (const auto& orig : spectrumWithoutErrors) {
        EXPECT_NE(std::find(spectrumWithErrors.begin(), spectrumWithErrors.end(), orig),
                  spectrumWithErrors.end());
    }
}

// Test 7: Sprawdzenie, że wynikowy spekturm jest posortowany leksykograficznie.
TEST_F(SpectrumGeneratorTest, GenerateSpectrumSortedOrder) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    std::string dna = "TACG"; // długość 4, k=2, bez błędów
    int k = 2;
    int errorCount = 0;
    // k-mery: "TA", "AC", "CG" – oczekiwany wynik posortowany: "AC", "CG", "TA"
    std::vector<std::string> expected = {"AC", "CG", "TA"};
    std::vector<std::string> spectrum = sg.generateSpectrum(dna, k, errorCount);
    EXPECT_EQ(spectrum, expected);
}

// Test 8: Sprawdzenie, że dla DNA z powtarzającymi się k-merami wynikowy spekturm zawiera tylko unikalne k-mery.
TEST_F(SpectrumGeneratorTest, GenerateSpectrumDuplicateHandling) {
    std::shared_ptr<std::mt19937> rng = std::make_shared<std::mt19937>(12345);
    SpectrumGenerator sg(rng);
    // DNA: "AAAAAA" (długość 6), k = 3 → wszystkie k-mery to "AAA"
    std::string dna = "AAAAAA";
    int k = 3;
    int errorCount = 0;
    std::vector<std::string> expected = {"AAA"};
    std::vector<std::string> spectrum = sg.generateSpectrum(dna, k, errorCount);
    EXPECT_EQ(spectrum, expected);
}