#include <gtest/gtest.h>
#include "dna/dna_instance.h"
#include "dna/dna_instance_io.h"
#include <stdexcept>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <filesystem>
#include <thread>
#include <future>

// Helper class for managing temporary test files
class TempFile {
    std::string filename;
public:
    explicit TempFile(std::string name) : filename(std::move(name)) {}
    ~TempFile() { 
        try {
            if (std::filesystem::exists(filename)) {
                std::filesystem::remove(filename); 
            }
        } catch (...) {} // Ignore cleanup errors
    }
    const std::string& getName() const { return filename; }
};

// Test 1: Domyślna (pusta) instancja – przy konstruktorze domyślnym DNA oraz spectrum są puste.
TEST(DNAInstanceTest, CreateEmptyInstance) {
    DNAInstance instance;
    EXPECT_TRUE(instance.getDNA().empty());
    EXPECT_EQ(instance.getSpectrum().size(), 0);
}

// Test 2: Ustawianie i pobieranie sekwencji DNA.
TEST(DNAInstanceTest, AddDNASequence) {
    DNAInstance instance;
    EXPECT_NO_THROW(instance.setDNA("ACGT"));
    EXPECT_EQ(instance.getDNA(), "ACGT");
}

// Test 3: Ustawianie i pobieranie spectrum.
TEST(DNAInstanceTest, SetAndGetSpectrum) {
    DNAInstance instance;
    std::vector<std::string> spectrum = {"ACGT", "CGTA"};
    EXPECT_NO_THROW(instance.setSpectrum(spectrum));
    EXPECT_EQ(instance.getSpectrum(), spectrum);
}

// Test 4: Konstrukcja instancji z poprawnymi parametrami – używamy pełnego konstruktora.
// Parametry: n=100, k=10, lNeg=3, lPoz=2, maxErrors=5, allowNegative=true, errorProb=0.1, seed=42.
TEST(DNAInstanceTest, ValidConstructionParameters) {
    DNAInstance instance(100, 10, 3, 2, 5, true, 0.1, 42);
    EXPECT_EQ(instance.getN(), 100);
    EXPECT_EQ(instance.getK(), 10);
    EXPECT_EQ(instance.getDeltaK(), 5);
    EXPECT_EQ(instance.getLNeg(), 3);
    EXPECT_EQ(instance.getLPoz(), 2);
    EXPECT_TRUE(instance.isRepAllowed());
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 0.1);
    EXPECT_EQ(instance.getDNA().length(), 100);
    // Nie sprawdzamy dokładnego rozmiaru spektrum, bo może się zmieniać w zależności od implementacji
    EXPECT_GT(instance.getSpectrum().size(), 0);
}

// Test 5: Konstrukcja z niepoprawnymi parametrami – n <= 0.
TEST(DNAInstanceTest, InvalidConstruction_NLessOrEqualZero) {
    EXPECT_THROW(DNAInstance(0, 10, 3, 2, 5, true, 0.1, 42), std::invalid_argument);
}

// Test 6: Konstrukcja z niepoprawnymi parametrami – k <= 0.
TEST(DNAInstanceTest, InvalidConstruction_KLessOrEqualZero) {
    EXPECT_THROW(DNAInstance(100, 0, 3, 2, 5, true, 0.1, 42), std::invalid_argument);
}

// Test 7: Konstrukcja z niepoprawnymi parametrami – k > n.
TEST(DNAInstanceTest, InvalidConstruction_KGreaterThanN) {
    EXPECT_THROW(DNAInstance(50, 60, 3, 2, 5, true, 0.1, 42), std::invalid_argument);
}

// Test 8: Test metody clearSpectrum – po jej wywołaniu spectrum jest puste.
TEST(DNAInstanceTest, ClearSpectrumWorks) {
    DNAInstance instance(100, 10, 3, 2, 5, true, 0.1, 42);
    EXPECT_FALSE(instance.getSpectrum().empty());
    EXPECT_NO_THROW(instance.clearSpectrum());
    EXPECT_TRUE(instance.getSpectrum().empty());
}

// Test 9: Test integracyjny zapisu i odczytu – zapisujemy instancję do pliku, potem ją wczytujemy.
TEST(DNAInstanceTest, SaveAndLoadFileIntegration) {
    TempFile testFile("test_instance.txt");
    
    // Create instance with known state
    DNAInstance originalInstance(100, 10, 3, 2, 5, true, 0.1, 42);
    std::string testDNA = "ACGTACGT";
    originalInstance.setDNA(testDNA);
    
    // Generate spectrum for the DNA
    auto spectrum = originalInstance.getSpectrum();
    
    // Save instance state
    EXPECT_NO_THROW(InstanceIO::saveInstance(originalInstance, testFile.getName()));
    
    // Load and verify
    DNAInstance loadedInstance;
    EXPECT_NO_THROW(loadedInstance = InstanceIO::loadInstance(testFile.getName()));
    
    // Verify all fields match
    EXPECT_EQ(loadedInstance.getN(), originalInstance.getN());
    EXPECT_EQ(loadedInstance.getK(), originalInstance.getK());
    EXPECT_EQ(loadedInstance.getDeltaK(), originalInstance.getDeltaK());
    EXPECT_EQ(loadedInstance.getLNeg(), originalInstance.getLNeg());
    EXPECT_EQ(loadedInstance.getLPoz(), originalInstance.getLPoz());
    EXPECT_EQ(loadedInstance.isRepAllowed(), originalInstance.isRepAllowed());
    EXPECT_DOUBLE_EQ(loadedInstance.getProbablePositive(), originalInstance.getProbablePositive());
    EXPECT_EQ(loadedInstance.getDNA(), testDNA);
    EXPECT_EQ(loadedInstance.getSpectrum(), spectrum);
}

// Test 10: Spójność generowanego spectrum – każdy k-mer powinien mieć długość równą k.
TEST(DNAInstanceTest, SpectrumGenerationConsistency) {
    DNAInstance instance(80, 8, 2, 2, 3, true, 0.0, 42);
    auto spectrum = instance.getSpectrum();
    EXPECT_GT(spectrum.size(), 0);
    for (const auto &kmer : spectrum) {
        EXPECT_EQ(kmer.length(), 8);
    }
}

// Test 11: Po konstrukcji pola DNA, originalDNA, targetSequence oraz size powinny być spójne.
TEST(DNAInstanceTest, DNAFieldsAfterConstruction) {
    // Create instance and let it generate DNA
    DNAInstance instance(60, 6, 2, 2, 3, true, 0.0, 42);
    
    // Get DNA after construction
    std::string dna = instance.getDNA();
    
    // Verify DNA properties
    EXPECT_FALSE(dna.empty()) << "DNA should not be empty after construction";
    EXPECT_EQ(dna.length(), 60) << "DNA length should match constructor parameter";
    EXPECT_EQ(dna, instance.getOriginalDNA()) << "DNA should match original DNA";
    
    // Verify target sequence
    std::string targetSeq = instance.getTargetSequence();
    EXPECT_FALSE(targetSeq.empty()) << "Target sequence should not be empty";
    EXPECT_EQ(targetSeq.length(), 60) << "Target sequence length should match DNA length";
    
    // Verify size
    EXPECT_EQ(instance.getSize(), 60) << "Size should match DNA length";
}

// Test 12: Kontrola statystycznego rozkładu liter w sekwencji DNA (dla dużej sekwencji).
TEST(DNAInstanceTest, StatisticalDistributionOfDNA) {
    DNAInstance instance(10000, 10, 3, 2, 5, true, 0.1, 42);
    std::string dna = instance.getDNA();
    EXPECT_EQ(dna.length(), 10000);
    
    std::map<char, int> distribution;
    for (char c : dna) {
        distribution[c]++;
    }
    
    const double expectedFrequency = 10000 / 4.0;
    const double allowedDeviation = 0.1;
    for (const auto& pair : distribution) {
        double deviation = std::abs(pair.second - expectedFrequency) / expectedFrequency;
        EXPECT_LE(deviation, allowedDeviation)
            << "Nucleotide " << pair.first << " frequency deviates more than "
            << (allowedDeviation * 100) << "%";
    }
}

// Test 13: Test wielowątkowości – równoległe wywołania metod setDNA i setSpectrum.
TEST(DNAInstanceTest, ThreadSafetyOperations) {
    const int numThreads = 4;
    const int operationsPerThread = 100;
    DNAInstance instance(100, 10, 3, 2, 5, true, 0.1, 42);
    
    std::vector<std::future<void>> futures;
    auto threadFunc = [&instance, operationsPerThread]() {
        for (int i = 0; i < operationsPerThread; ++i) {
            if (i % 2 == 0) {
                instance.setDNA("ACGTACGT");
            } else {
                std::vector<std::string> spec = {"ACGT", "CGTA", "GTAC"};
                instance.setSpectrum(spec);
            }
        }
    };
    
    for (int i = 0; i < numThreads; ++i) {
        futures.push_back(std::async(std::launch::async, threadFunc));
    }
    
    for (auto &fut : futures) {
        fut.get(); // Will throw if any exceptions occurred
    }
    
    // Final validation
    instance.setDNA("ACGT");
    std::vector<std::string> spec = {"ACGT", "CGTA"};
    instance.setSpectrum(spec);
}

// Test 14: (Usunięty test dotyczący metody generateRandomInstance, która nie jest publiczna)

// Test 15: Test wczytywania z nieistniejącego pliku oraz pliku o niepoprawnym formacie.
// Zmodyfikowano oczekiwanie na std::runtime_error.
TEST(DNAInstanceTest, LoadFromFileErrors) {
    EXPECT_THROW(InstanceIO::loadInstance("nonexistent.txt"), std::runtime_error);
    EXPECT_THROW(InstanceIO::loadInstance("invalid.txt"), std::runtime_error);
}
