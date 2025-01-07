//
// Created by konrad_guest on 29/12/2024.
//

#ifndef NAIVE_RECONSTRUCTION_H
#define NAIVE_RECONSTRUCTION_H
#include <string>
#include <vector>
#include "generator/dna_generator.h"


/**
 * Enum z trzema metodami naiwnej rekonstrukcji.
 * Możesz nazwać je inaczej lub dodać kolejne.
 */
enum class NaiveReconstructionMethod {
    METHOD_A,  // np. sortowanie alfabetyczne i sklejanie
    METHOD_B,  // np. łączenie sufiksów i prefiksów
    METHOD_C   // np. inny prosty heurystyczny pomysł
};

/**
 * Klasa zbiorcza, oferująca 3 różne (naiwne) sposoby
 * na zrekonstruowanie DNA z danego spektrum.
 */
class NaiveReconstructor {
public:
    /**
     * Metoda główna: na podstawie "method" wywołuje jedną z trzech
     * implementacji w tej klasie.
     */
    std::string reconstructDNA(const DNAInstance &instance,
                               NaiveReconstructionMethod method);

private:
    // Każda z poniższych metod realizuje inny (naiwny) pomysł.

    /**
     * Metoda A: np. sortowanie k-merów alfabetycznie i proste ich sklejanie.
     */
    std::string reconstructDNA_A(const DNAInstance &instance);

    /**
     * Metoda B: np. łączenie k-merów za pomocą minimalnego dopasowania
     * sufiks -> prefiks (w bardzo naiwny sposób).
     */
    std::string reconstructDNA_B(const DNAInstance &instance);

    /**
     * Metoda C: np. bierzemy losowo k-mery i sklejamy, doklejając po (k-1) znaków,
     * albo inne dowolne, proste rozwiązanie.
     */
    std::string reconstructDNA_C(const DNAInstance &instance);
};


/**
 * Funkcja obliczająca odległość Levenshteina między dwoma łańcuchami znaków.
 * Zwraca najmniejszą liczbę operacji wstawiania, usuwania lub zamiany
 * potrzebnych do przekształcenia `s1` w `s2`.
 */
int levenshteinDistance(const std::string &s1, const std::string &s2);



#endif //NAIVE_RECONSTRUCTION_H
