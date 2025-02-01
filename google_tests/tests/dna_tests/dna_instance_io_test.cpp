//
// Created by konrad_guest on 01/02/2025.
//
#include <gtest/gtest.h>
#include "dna/dna_instance_io.h"
#include "dna/dna_instance.h"
#include "utils/logging.h"
#include "../base_test.h"  // Bazowy plik testowy (zapewnia np. SetUp/TearDown logera)
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

// Klasa pomocnicza – dziedziczy po BaseTest, która inicjalizuje logera
class DNAInstanceIOTest : public BaseTest {
protected:
    void SetUp() override {
        Logger::initialize("dna_instance_io_test.log");
        Logger::setLogLevel(LogLevel::INFO);
        BaseTest::SetUp();
    }

    void TearDown() override {
        BaseTest::TearDown();
        Logger::cleanup();
    }

    // Funkcja pomocnicza do utworzenia przykładowej instancji DNA
    DNAInstance createTestInstance() {
        DNAInstance instance;
        // Ustawiamy parametry – przykładowe wartości
        instance.setN(100);
        instance.setK(10);
        instance.setDeltaK(2);
        instance.setLNeg(5);
        instance.setLPoz(5);
        instance.setRepAllowed(true);
        instance.setProbablePositive(0.0);
        instance.setSize(100);
        instance.setDNA("ACGTACGTAC");
        instance.setOriginalDNA("ACGTACGTAC");
        instance.setTargetSequence("ACGTACGTAC");
        std::vector<std::string> spectrum = {"ACGTACGTAC", "CGTACGTACG", "GTACGTACGT", "TACGTACGTA"};
        instance.setSpectrum(spectrum);
        return instance;
    }
};

// Test 1: Zapis i odczyt poprawnej instancji – po zapisie instancja jest taka sama, jak przed zapisem.
TEST_F(DNAInstanceIOTest, SaveAndLoadInstanceSuccess) {
    DNAInstance instance = createTestInstance();
    std::string filename = "test_instance_io.txt";

    // Zapis instancji – nie powinno być wyjątku i metoda powinna zwrócić true.
    EXPECT_NO_THROW({
        bool saved = InstanceIO::saveInstance(instance, filename);
        EXPECT_TRUE(saved);
    });

    // Odczytujemy instancję
    DNAInstance loaded;
    EXPECT_NO_THROW({
        loaded = InstanceIO::loadInstance(filename);
    });

    // Weryfikacja: porównujemy wszystkie istotne pola
    EXPECT_EQ(loaded.getN(), instance.getN());
    EXPECT_EQ(loaded.getK(), instance.getK());
    EXPECT_EQ(loaded.getDeltaK(), instance.getDeltaK());
    EXPECT_EQ(loaded.getLNeg(), instance.getLNeg());
    EXPECT_EQ(loaded.getLPoz(), instance.getLPoz());
    EXPECT_EQ(loaded.isRepAllowed(), instance.isRepAllowed());
    EXPECT_EQ(loaded.getProbablePositive(), instance.getProbablePositive());
    EXPECT_EQ(loaded.getSize(), instance.getSize());
    EXPECT_EQ(loaded.getDNA(), instance.getDNA());
    EXPECT_EQ(loaded.getOriginalDNA(), instance.getOriginalDNA());
    EXPECT_EQ(loaded.getTargetSequence(), instance.getTargetSequence());
    EXPECT_EQ(loaded.getSpectrum(), instance.getSpectrum());

    // Usuwamy plik testowy
    std::filesystem::remove(filename);
}

// Test 2: Próba zapisu przy pustej nazwie pliku – metoda powinna rzucić wyjątek.
TEST_F(DNAInstanceIOTest, SaveInstanceInvalidFilename) {
    DNAInstance instance = createTestInstance();
    EXPECT_THROW({
        InstanceIO::saveInstance(instance, "");
    }, std::runtime_error);
}

// Test 3: Próba odczytu z nieistniejącego pliku – metoda powinna rzucić wyjątek.
TEST_F(DNAInstanceIOTest, LoadInstanceNonexistentFile) {
    EXPECT_THROW({
        InstanceIO::loadInstance("nonexistent_file.txt");
    }, std::runtime_error);
}

// Test 4: Odczyt z uszkodzonego pliku – najpierw zapisujemy poprawną instancję,
// następnie modyfikujemy plik (np. skracamy zawartość), a potem oczekujemy wyjątku.
TEST_F(DNAInstanceIOTest, LoadInstanceCorruptedFile) {
    DNAInstance instance = createTestInstance();
    std::string filename = "test_corrupted_instance.txt";

    // Zapisujemy instancję
    EXPECT_NO_THROW({
        bool saved = InstanceIO::saveInstance(instance, filename);
        EXPECT_TRUE(saved);
    });

    // Korumpujemy plik – otwieramy i zastępujemy całą zawartość
    std::ofstream file(filename, std::ios::out | std::ios::trunc);
    ASSERT_TRUE(file.is_open());
    file << "corrupted data";
    file.close();

    // Próba odczytu powinna rzucić wyjątkiem
    EXPECT_THROW({
        InstanceIO::loadInstance(filename);
    }, std::runtime_error);

    std::filesystem::remove(filename);
}

// Test 5: Sprawdzenie formatu zapisu – odczytujemy plik tekstowy i sprawdzamy, czy pierwszy wiersz
// zawiera wartość odpowiadającą getN() instancji.
TEST_F(DNAInstanceIOTest, SaveInstanceFileContentFormat) {
    DNAInstance instance = createTestInstance();
    std::string filename = "test_instance_format.txt";

    EXPECT_NO_THROW({
        bool saved = InstanceIO::saveInstance(instance, filename);
        EXPECT_TRUE(saved);
    });

    // Odczytujemy zawartość pliku jako tekst
    std::ifstream file(filename);
    ASSERT_TRUE(file.is_open());
    std::stringstream ss;
    ss << file.rdbuf();
    std::string content = ss.str();
    file.close();

    // Pobieramy pierwszy wiersz i porównujemy z wartością getN()
    std::istringstream iss(content);
    std::string firstLine;
    std::getline(iss, firstLine);
    EXPECT_EQ(firstLine, std::to_string(instance.getN()));

    std::filesystem::remove(filename);
}

// Test 6: Próba odczytu pliku niekompletnego – plik, w którym brakuje części danych (np. spektrum)
// powinien spowodować rzucenie wyjątku.
TEST_F(DNAInstanceIOTest, LoadInstanceIncompleteFile) {
    std::string filename = "test_incomplete_instance.txt";
    std::ofstream file(filename);
    ASSERT_TRUE(file.is_open());
    // Zapisujemy tylko podstawowe parametry i DNA, pomijając spektrum
    file << "100\n";   // n
    file << "10\n";    // k
    file << "2\n";     // deltaK
    file << "5\n";     // lNeg
    file << "5\n";     // lPoz
    file << "1\n";     // repAllowed
    file << "0.0\n";   // probablePositive
    file << "100\n";   // size
    file << "ACGTACGTAC\n";  // DNA
    file << "ACGTACGTAC\n";  // originalDNA
    file << "ACGTACGTAC\n";  // targetSequence
    // Brak zapisu spekturm (nie ma liczby k-merów ani ich zawartości)
    file.close();

    EXPECT_THROW({
        InstanceIO::loadInstance(filename);
    }, std::runtime_error);

    std::filesystem::remove(filename);
}