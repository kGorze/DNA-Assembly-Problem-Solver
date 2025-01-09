//
// Created by konrad_guest on 07/01/2025.
//

#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include <string>
#include <vector>
#include "generator/dna_generator.h"

/**
 * Typ wyliczeniowy – wybór reprezentacji osobnika.
 */
enum class RepresentationType {
    DIRECT,      // A: Łańcuch DNA długości n
    PERMUTATION, // B: Permutacja k-merów
    GRAPH_PATH   // C: Ścieżka w grafie nakładających się k-merów
};

/**
 * Interfejs wspólny dla reprezentacji osobnika w algorytmie genetycznym.
 * Pozwoli nam na jednolite wywołania w GeneticAlgorithm i Fitness.
 */
class IRepresentation {
public:
    virtual ~IRepresentation() = default;

    /**
     * Inicjalizacja populacji – zależnie od typu reprezentacji
     */
    virtual std::vector<void*> initializePopulation(int popSize, const DNAInstance &instance) = 0;

    /**
     * Krótkie wytłumaczenie:
     * - "Osobnik" może być w różnej postaci (np. std::string, std::vector<int>, 
     *   a w przypadku grafu – pewna struktura ścieżki).
     * - Trzymamy go jako `void*` w populacji, by nie ograniczać się do jednego typu.
     * - W praktyce lepiej użyć std::variant lub polimorfizmu, ale tutaj pokazujemy ideę.
     */

    /**
     * Ewentualna funkcja dekodująca osobnika do postaci "łańcucha DNA" – 
     * potrzebna np. w fitness, gdy chcemy oceniać liczbę dopasowanych k-merów
     * albo odległość Levenshteina do oryginału.
     */
    virtual std::string decodeToDNA(void* individual, const DNAInstance &instance) = 0;
};

/**
 * Klasa A: Bezpośredni łańcuch DNA – "DirectDNARepresentation"
 */
class DirectDNARepresentation : public IRepresentation {
public:
  DirectDNARepresentation(int fallbackN = 500, double randomRatio = 1.0)
    : fallbackN(fallbackN), randomRatio(randomRatio) {}
  
    std::vector<void*> initializePopulation(int popSize, const DNAInstance &instance) override;

 
    std::string decodeToDNA(void* individual, const DNAInstance &instance) override;

private:
  int fallbackN;
  double randomRatio;
};

/**
 * Klasa B: Permutacja k-merów – "PermutationRepresentation"
 */
class PermutationRepresentation : public IRepresentation {
public:
    std::vector<void*> initializePopulation(int popSize, const DNAInstance &instance) override;
    std::string decodeToDNA(void* individual, const DNAInstance &instance) override;
};

/**
 * Klasa C: Ścieżka w grafie nakładania k-merów – "GraphPathRepresentation"
 * (uwaga: przykład w dużym uproszczeniu!)
 */
class GraphPathRepresentation : public IRepresentation {
public:
    std::vector<void*> initializePopulation(int popSize, const DNAInstance &instance) override;
    std::string decodeToDNA(void* individual, const DNAInstance &instance) override;
};

#endif //REPRESENTATION_H
