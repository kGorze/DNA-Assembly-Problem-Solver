//
// Created by konrad_guest on 01/02/2025.
//
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
#include <string>
#include <memory>

// Uwzględnij odpowiednie ścieżki – dostosuj do swojej struktury projektu:
#include "dna/dna_instance_builder.h"
#include "generator/dna_generator.h"
#include "dna/error_introduction.h"
#include "utils/random.h"

// ----- Dummy strategie błędu -----

// Dummy strategy – usuwa pierwszy k-mer ze spectrum
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

// Strategy, która zawsze rzuca wyjątek
class FailingErrorStrategy : public IErrorIntroductionStrategy {
public:
    void introduceErrors(DNAInstance& instance) override {
        (void)instance; // unikamy ostrzeżenia o nieużywanym parametrze
        throw std::runtime_error("Intentional error in error strategy");
    }
};

// Funkcja pomocnicza – tworzy poprawny generator DNA
std::shared_ptr<DNAGenerator> createValidGenerator() {
    // Konstrukcja DNAGenerator wymaga przekazania unikalnego obiektu Random
    return std::make_shared<DNAGenerator>(std::make_unique<Random>());
}

// ==================== TESTY DNAInstanceBuilder ====================

// 1. Inicjalizacja buildera z poprawnym generatorem nie rzuca wyjątku.
TEST(DNAInstanceBuilderTest, InitializationWithValidGenerator) {
    auto generator = createValidGenerator();
    EXPECT_NO_THROW({
        DNAInstanceBuilder builder(generator);
    });
}

// 2. Inicjalizacja buildera z nullptr jako generator – powinno rzucić wyjątek.
TEST(DNAInstanceBuilderTest, InitializationWithNullGenerator) {
    EXPECT_THROW({
        DNAInstanceBuilder builder(nullptr);
    }, std::invalid_argument);
}

// 3. Ustawianie N poprzez builder.
TEST(DNAInstanceBuilderTest, SetNValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(100);
    EXPECT_EQ(builder.getInstance().getN(), 100);
}

// 4. Ustawianie K poprzez builder.
TEST(DNAInstanceBuilderTest, SetKValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setK(10);
    EXPECT_EQ(builder.getInstance().getK(), 10);
}

// 5. Ustawianie deltaK.
TEST(DNAInstanceBuilderTest, SetDeltaKValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setDeltaK(2);
    EXPECT_EQ(builder.getInstance().getDeltaK(), 2);
}

// 6. Ustawianie LNeg.
TEST(DNAInstanceBuilderTest, SetLNegValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setLNeg(5);
    EXPECT_EQ(builder.getInstance().getLNeg(), 5);
}

// 7. Ustawianie LPoz.
TEST(DNAInstanceBuilderTest, SetLPozValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setLPoz(3);
    EXPECT_EQ(builder.getInstance().getLPoz(), 3);
}

// 8. Ustawianie flagi repAllowed.
TEST(DNAInstanceBuilderTest, SetRepAllowedValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setRepAllowed(false);
    EXPECT_FALSE(builder.getInstance().isRepAllowed());
}

// 9. Ustawianie probablePositive.
TEST(DNAInstanceBuilderTest, SetProbablePositiveValid) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setProbablePositive(0.25);
    EXPECT_DOUBLE_EQ(builder.getInstance().getProbablePositive(), 0.25);
}

// 10. Test buildDNA – po buildDNA sekwencja DNA nie jest pusta i ma długość równą N.
TEST(DNAInstanceBuilderTest, BuildDNAGeneratesDNA) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(50).setK(5).setDeltaK(1);
    builder.buildDNA();
    std::string dna = builder.getInstance().getDNA();
    EXPECT_FALSE(dna.empty());
    EXPECT_EQ(dna.length(), 50);
}

// 11. Test buildSpectrum – po buildSpectrum spectrum nie jest puste, a każdy k-mer ma długość równą K.
TEST(DNAInstanceBuilderTest, BuildSpectrumGeneratesSpectrum) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(50).setK(5).setDeltaK(1).buildDNA().buildSpectrum();
    auto spectrum = builder.getInstance().getSpectrum();
    EXPECT_FALSE(spectrum.empty());
    for (const auto &kmer : spectrum) {
        EXPECT_EQ(kmer.length(), 5);
    }
}

// 12. Pełny łańcuch: ustawienie parametrów i budowanie instancji bez błędów.
TEST(DNAInstanceBuilderTest, FullChainTestWithoutErrors) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80)
           .setK(7)
           .setDeltaK(1)
           .setLNeg(3)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.0)
           .buildDNA()
           .buildSpectrum();
    DNAInstance instance = builder.getInstance();
    EXPECT_EQ(instance.getN(), 80);
    EXPECT_EQ(instance.getK(), 7);
    EXPECT_GE(instance.getDeltaK(), 0);
    EXPECT_EQ(instance.getLNeg(), 3);
    EXPECT_EQ(instance.getLPoz(), 2);
    EXPECT_TRUE(instance.isRepAllowed());
    EXPECT_DOUBLE_EQ(instance.getProbablePositive(), 0.0);
    EXPECT_FALSE(instance.getDNA().empty());
    EXPECT_FALSE(instance.getSpectrum().empty());
}

// 13. Test applyError z nullptr – metoda nie powinna rzucać wyjątku.
TEST(DNAInstanceBuilderTest, ApplyErrorWithNullStrategy) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80).setK(7).setDeltaK(1).buildDNA().buildSpectrum();
    EXPECT_NO_THROW({
        builder.applyError(nullptr);
    });
}

// 14. Test applyError z DummyErrorStrategy – spectrum powinno się zmniejszyć o jeden element.
TEST(DNAInstanceBuilderTest, ApplyErrorWithDummyStrategy) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80).setK(7).setDeltaK(1).buildDNA().buildSpectrum();
    auto before = builder.getInstance().getSpectrum();
    DummyErrorStrategy strategy;
    builder.applyError(&strategy);
    auto after = builder.getInstance().getSpectrum();
    EXPECT_EQ(after.size(), before.size() - 1);
}

// 15. Test wielokrotnych aktualizacji parametrów – ostatnia wartość ma być zachowana.
TEST(DNAInstanceBuilderTest, MultipleParameterUpdates) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(100).setK(10).setDeltaK(1);
    DNAInstance inst1 = builder.getInstance();
    EXPECT_EQ(inst1.getN(), 100);
    builder.setN(150);  // aktualizacja
    DNAInstance inst2 = builder.getInstance();
    EXPECT_EQ(inst2.getN(), 150);
}

// 16. Test jawnego ustawiania sekwencji DNA przez setDNA.
TEST(DNAInstanceBuilderTest, SetDNAExplicitly) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    std::string customDNA = "ACGTACGTACGT";
    builder.setDNA(customDNA);
    EXPECT_EQ(builder.getInstance().getDNA(), customDNA);
}

// 17. Test jawnego ustawiania spectrum przez setSpectrum.
TEST(DNAInstanceBuilderTest, SetSpectrumExplicitly) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    std::vector<std::string> customSpectrum = {"ACG", "CGT", "GTA"};
    builder.setSpectrum(customSpectrum);
    EXPECT_EQ(builder.getInstance().getSpectrum(), customSpectrum);
}

// 18. Test – jeżeli nie ustawimy wymaganych parametrów (N lub K) – buildDNA powinno rzucić wyjątkiem.
TEST(DNAInstanceBuilderTest, BuildDNAWithoutMandatoryParameters) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    // Nie ustawiamy N ani K – oczekujemy wyjątku przy wywołaniu buildDNA.
    EXPECT_THROW({
        builder.buildDNA();
    }, std::invalid_argument);
}

// 19. Test integracyjny: pełny łańcuch (buildDNA, buildSpectrum, applyError, build) daje spójną instancję.
TEST(DNAInstanceBuilderTest, FullChainBuildAndApplyErrorIntegration) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(90)
           .setK(7)
           .setDeltaK(1)
           .setLNeg(2)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.1)
           .buildDNA()
           .buildSpectrum()
           .applyError(nullptr); // brak dodatkowych błędów
    EXPECT_NO_THROW({
        DNAInstance finalInstance = builder.build();
        EXPECT_FALSE(finalInstance.getDNA().empty());
        EXPECT_FALSE(finalInstance.getSpectrum().empty());
    });
}

// 20. Test walidacji stanu w kontekście "buildDNA" – przed generowaniem spectrum metoda buildDNA powinna działać.
TEST(DNAInstanceBuilderTest, ValidateStateContextBuildDNA) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80).setK(7).setDeltaK(1);
    EXPECT_NO_THROW({
        builder.buildDNA();
    });
}

// 21. Test walidacji stanu w kontekście "buildSpectrum" – jeżeli spectrum jest puste, buildSpectrum powinno rzucić wyjątek.
// Aby wymusić pusty spectrum, używamy const_cast.
TEST(DNAInstanceBuilderTest, ValidateStateContextBuildSpectrumEmptySpectrum) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80).setK(7).setDeltaK(1).buildDNA();
    // Uzyskaj nie-const referencję poprzez const_cast
    DNAInstance &inst = const_cast<DNAInstance&>(builder.getInstance());
    inst.clearSpectrum();
    EXPECT_THROW({
        builder.buildSpectrum();
    }, std::runtime_error);
}

// 22. Test – po buildSpectrum wszystkie k-mery mają długość równą K.
TEST(DNAInstanceBuilderTest, KmerLengthConsistencyAfterBuildSpectrum) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(100).setK(5).setDeltaK(1).buildDNA().buildSpectrum();
    auto spectrum = builder.getInstance().getSpectrum();
    for (const auto &kmer : spectrum) {
        EXPECT_EQ(kmer.length(), 5);
    }
}

// 23. Test – po buildDNA pola DNA, originalDNA, targetSequence oraz size są poprawnie ustawione.
TEST(DNAInstanceBuilderTest, DNAFieldsAfterBuildDNA) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(60).setK(6).setDeltaK(1).buildDNA();
    DNAInstance instance = builder.getInstance();
    EXPECT_FALSE(instance.getDNA().empty());
    EXPECT_EQ(instance.getDNA(), instance.getOriginalDNA());
    EXPECT_EQ(instance.getDNA(), instance.getTargetSequence());
    EXPECT_EQ(instance.getSize(), instance.getDNA().length());
}

// 24. Test – finalna instancja (build) musi mieć ustawione wszystkie wymagane pola.
TEST(DNAInstanceBuilderTest, FinalInstanceIntegrity) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80)
           .setK(7)
           .setDeltaK(1)
           .setLNeg(2)
           .setLPoz(2)
           .setRepAllowed(true)
           .setProbablePositive(0.0)
           .buildDNA()
           .buildSpectrum();
    EXPECT_NO_THROW({
        DNAInstance finalInstance = builder.build();
        EXPECT_FALSE(finalInstance.getDNA().empty());
        EXPECT_FALSE(finalInstance.getSpectrum().empty());
    });
}

// 25. Test – applyError z FailingErrorStrategy powinno rzucić wyjątek.
TEST(DNAInstanceBuilderTest, ApplyErrorExceptionHandling) {
    auto generator = createValidGenerator();
    DNAInstanceBuilder builder(generator);
    builder.setN(80).setK(7).setDeltaK(1).buildDNA().buildSpectrum();
    FailingErrorStrategy failingStrategy;
    EXPECT_THROW({
        builder.applyError(&failingStrategy);
    }, std::runtime_error);
}
