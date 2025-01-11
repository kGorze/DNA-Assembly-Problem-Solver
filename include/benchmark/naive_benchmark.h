//
// Created by konrad_guest on 29/12/2024.
// SMART

#ifndef NAIVE_BENCHMARK_H
#define NAIVE_BENCHMARK_H
#include "naive/naive_reconstruction.h"
#include <chrono>

class NaiveBenchmark {
public:
    /**
     * Uruchamia testy (benchmark) dla wszystkich trzech metod
     * z klasy NaiveReconstructor.
     *
     * Możesz tu zdefiniować dowolną strukturę zwracaną
     * (np. JSON, struct z wynikami itp.).
     */
    void runBenchmark(const DNAInstance &instance);
};

#endif //NAIVE_BENCHMARK_H
