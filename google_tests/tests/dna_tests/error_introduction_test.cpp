//
// Created by konrad_guest on 01/02/2025.
//
#include <gtest/gtest.h>
#include "dna/error_introduction.h"
#include "dna/dna_instance.h"
#include "utils/logging.h"
#include "../base_test.h"  // Zakładamy, że BaseTest zawiera podstawową konfigurację środowiska testowego
#include <filesystem>
#include <stdexcept>
#include <vector>
#include <string>
#include <algorithm>

// Funkcja pomocnicza tworząca przykładową, poprawną instancję DNA
// Zakładamy, że DNAInstance posiada odpowiednie metody set/get.
DNAInstance createValidDNAInstance() {
    DNAInstance instance;
    // Ustawiamy przykładowe parametry:
    instance.setN(100);
    instance.setK(10);
    instance.setDeltaK(2);
    instance.setLNeg(5);   // liczba błędów negatywnych
    instance.setLPoz(3);   // liczba błędów pozytywnych
    instance.setRepAllowed(true);
    instance.setProbablePositive(0);
    instance.setSize(100);

    // Ustawiamy przykładową sekwencję DNA – nie musi być realistyczna, ale niepusta.
    instance.setDNA("ACGTACGTAC");
    instance.setOriginalDNA("ACGTACGTAC");
    instance.setTargetSequence("ACGTACGTAC");

    // Ustawiamy przykładowy spekturm – zakładamy, że mamy kilka k-merów
    std::vector<std::string> spectrum = {
        "ACGTACGTAC", "CGTACGTACG", "GTACGTACGT", "TACGTACGTA", "ACGTACGTAC"
    };
    instance.setSpectrum(spectrum);

    return instance;
}

// Funkcja pomocnicza tworząca instancję o pustym DNA (niepoprawna)
DNAInstance createEmptyDNADNAInstance() {
    DNAInstance instance = createValidDNAInstance();
    instance.setDNA("");
    return instance;
}

// Funkcja pomocnicza tworząca instancję o pustym spekturm-ie (niepoprawna)
DNAInstance createEmptySpectrumInstance() {
    DNAInstance instance = createValidDNAInstance();
    std::vector<std::string> emptySpec;
    instance.setSpectrum(emptySpec);
    return instance;
}

// Klasa testowa dziedzicząca po BaseTest – inicjalizacja logera
class ErrorIntroductionTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("error_introduction_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();
    }
    void TearDown() override {
        BaseTest::TearDown();
        Logger::cleanup();
    }
};

// ----- NegativeErrorIntroducer Tests -----

// Test 1: Jeśli DNA jest puste, metoda validateInstance powinna zwrócić false, a introduceErrors rzucić wyjątek.
TEST_F(ErrorIntroductionTest, NegativeError_InvalidDNA) {
    DNAInstance instance;  // Tworzymy nową instancję bez ustawiania DNA
    NegativeErrorIntroducer negErr(5);
    EXPECT_THROW({
        negErr.introduceErrors(instance);
    }, std::runtime_error);
}

// Test 2: Jeśli spekturm jest pusty, metoda validateInstance powinna zwrócić false, a introduceErrors rzucić wyjątek.
TEST_F(ErrorIntroductionTest, NegativeError_InvalidSpectrum) {
    DNAInstance instance = createEmptySpectrumInstance();
    NegativeErrorIntroducer negErr(5);
    EXPECT_THROW({
        negErr.introduceErrors(instance);
    }, std::runtime_error);
}

// Test 3: Jeśli wartość lNeg <= 0, introduceErrors nie powinno zmieniać spekturm-u.
TEST_F(ErrorIntroductionTest, NegativeError_NoErrorsWhenLNegNonPositive) {
    DNAInstance instance = createValidDNAInstance();
    // Ustawiamy lNeg na 0
    instance.setLNeg(0);
    std::vector<std::string> originalSpectrum = instance.getSpectrum();

    NegativeErrorIntroducer negErr(0);
    // Nie powinno rzucić wyjątku, a spekturm pozostanie taki sam.
    EXPECT_NO_THROW({
        negErr.introduceErrors(instance);
    });
    EXPECT_EQ(instance.getSpectrum(), originalSpectrum);
}

// Test 4: Poprawne wprowadzenie błędów – po introduceErrors rozmiar spekturm-u zmniejsza się (ale nie jest pusty).
TEST_F(ErrorIntroductionTest, NegativeError_IntroduceErrors) {
    DNAInstance instance = createValidDNAInstance();
    int originalSize = instance.getSpectrum().size();
    // Upewniamy się, że lNeg jest mniejszy niż rozmiar spektrum
    instance.setLNeg(2);  // Ustawiamy na 2, bo wiemy że spektrum ma 5 elementów
    NegativeErrorIntroducer negErr(2);
    EXPECT_NO_THROW({
        negErr.introduceErrors(instance);
    });
    int newSize = instance.getSpectrum().size();
    // Spodziewamy się, że usunięto 2 elementy
    EXPECT_EQ(newSize, originalSize - 2);
    EXPECT_GT(newSize, 0);
}

// ----- PositiveErrorIntroducer Tests -----

// Test 5: Jeśli DNA jest puste, metoda validateInstance powinna zwrócić false, a introduceErrors rzucić wyjątek.
TEST_F(ErrorIntroductionTest, PositiveError_InvalidDNA) {
    DNAInstance instance;  // Tworzymy nową instancję bez ustawiania DNA
    PositiveErrorIntroducer posErr(3, 10);
    EXPECT_THROW({
        posErr.introduceErrors(instance);
    }, std::runtime_error);
}

// Test 6: Jeśli spekturm jest pusty, metoda validateInstance powinna zwrócić false, a introduceErrors rzucić wyjątek.
TEST_F(ErrorIntroductionTest, PositiveError_InvalidSpectrum) {
    DNAInstance instance = createEmptySpectrumInstance();
    PositiveErrorIntroducer posErr(3, instance.getK());
    EXPECT_THROW({
        posErr.introduceErrors(instance);
    }, std::runtime_error);
}

// Test 7: Jeśli wartość lPoz <= 0, introduceErrors nie powinno zmieniać spekturm-u.
TEST_F(ErrorIntroductionTest, PositiveError_NoErrorsWhenLPozNonPositive) {
    DNAInstance instance = createValidDNAInstance();
    // Ustawiamy LPoz na 0
    instance.setLPoz(0);
    std::vector<std::string> originalSpectrum = instance.getSpectrum();

    PositiveErrorIntroducer posErr(0, instance.getK());
    EXPECT_NO_THROW({
        posErr.introduceErrors(instance);
    });
    EXPECT_EQ(instance.getSpectrum(), originalSpectrum);
}

// Test 8: Poprawne wprowadzenie błędów pozytywnych – po introduceErrors rozmiar spekturm-u zwiększa się o wartość LPoz.
TEST_F(ErrorIntroductionTest, PositiveError_IntroduceErrors) {
    DNAInstance instance = createValidDNAInstance();
    int originalSize = instance.getSpectrum().size();
    // Ustaw LPoz na 3 (liczba błędów pozytywnych do wprowadzenia)
    instance.setLPoz(3);
    PositiveErrorIntroducer posErr(3, instance.getK());
    EXPECT_NO_THROW({
        posErr.introduceErrors(instance);
    });
    int newSize = instance.getSpectrum().size();
    // Oczekujemy, że do spekturm-u zostaną dodane 3 nowe unikalne k-mery.
    EXPECT_EQ(newSize, originalSize + 3);
    // Dodatkowo, wszystkie nowe k-mery powinny mieć długość równą getK()
    for (const auto& kmer : instance.getSpectrum()) {
        EXPECT_EQ(kmer.length(), static_cast<size_t>(instance.getK()));
    }
}

// Test 9: Test metody generateRandomKmer – sprawdzamy, czy zwraca ciąg o zadanej długości i zawierający tylko "ACGT".
TEST_F(ErrorIntroductionTest, PositiveError_GenerateRandomKmer) {
    PositiveErrorIntroducer posErr(3, 7);
    int length = 7;
    std::string kmer = posErr.generateRandomKmer(length);
    EXPECT_EQ(kmer.length(), static_cast<size_t>(length));
    std::string validChars = "ACGT";
    for (char c : kmer) {
        EXPECT_NE(validChars.find(c), std::string::npos)
            << "Wygenerowano nieprawidłowy znak: " << c;
    }
}

// Test 10: Test fabryki ErrorIntroducerFactory – sprawdzamy, czy zwracane obiekty są zgodne z interfejsem IErrorIntroductionStrategy.
TEST_F(ErrorIntroductionTest, FactoryCreatesErrorIntroducers) {
    auto negIntroducer = ErrorIntroducerFactory::createNegativeErrorIntroducer(4);
    auto posIntroducer = ErrorIntroducerFactory::createPositiveErrorIntroducer(4, 10);
    // Sprawdzamy, czy wskaźniki nie są puste (możemy także dynamic_cast, jeśli interfejs jest polimorficzny)
    EXPECT_NE(negIntroducer, nullptr);
    EXPECT_NE(posIntroducer, nullptr);
}

// Test 11: Sprawdzenie sytuacji, w której po usunięciu k-merów przez NegativeErrorIntroducer spekturm staje się pusty – powinna zostać rzucona odpowiednia sytuacja.
TEST_F(ErrorIntroductionTest, NegativeError_ResultingEmptySpectrumThrows) {
    DNAInstance instance = createValidDNAInstance();
    // Ustawiamy lNeg równy rozmiarowi spektrum
    std::vector<std::string> spectrum = {"ACGTACGTAC", "CGTACGTACG", "GTACGTACGT"};
    instance.setSpectrum(spectrum);
    instance.setLNeg(3);
    
    NegativeErrorIntroducer negErr(3);
    EXPECT_THROW({
        negErr.introduceErrors(instance);
    }, std::runtime_error);
}