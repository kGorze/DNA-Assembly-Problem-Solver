//
// Created by konrad_guest on 29/12/2024.
// SMART

#include "benchmark/naive_benchmark.h"

void NaiveBenchmark::runBenchmark(const DNAInstance &instance)
{
    NaiveReconstructor reconstructor;

    // Oryginalne DNA (z pliku lub pamięci)
    std::string original = instance.getDNA();

    // Lista metod do przetestowania
    NaiveReconstructionMethod methods[3] = {
        NaiveReconstructionMethod::METHOD_A,
        NaiveReconstructionMethod::METHOD_B,
        NaiveReconstructionMethod::METHOD_C
    };

    // Nazwy do wypisania
    const char* methodNames[3] = {"METHOD_A", "METHOD_B", "METHOD_C"};

    for(int i = 0; i < 3; ++i) {
        // Pomiar czasu
        auto t1 = std::chrono::high_resolution_clock::now();
        std::string reconstructed = reconstructor.reconstructDNA(instance, methods[i]);
        auto t2 = std::chrono::high_resolution_clock::now();

        // Oblicz czas
        double elapsed = std::chrono::duration<double, std::milli>(t2 - t1).count();

        // Oblicz odległość Levenshteina
        int dist = levenshteinDistance(original, reconstructed);

        // Wypisz wyniki – w realnym projekcie można zapisać do pliku/tabeli CSV
        std::cout << "[NaiveBenchmark] " << methodNames[i] << ": "
                  << "Czas = " << elapsed << " ms, "
                  << "Levenshtein = " << dist << "\n";
    }
}
