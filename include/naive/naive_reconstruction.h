//
// Created by konrad_guest on 29/12/2024.
// SMART

#ifndef NAIVE_RECONSTRUCTION_H
#define NAIVE_RECONSTRUCTION_H
#include <string>
#include <vector>
#include "generator/dna_generator.h"
#include "utils/utility_functions.h"

/**
 * Enum z trzema metodami naiwnej rekonstrukcji.
 */
enum class NaiveReconstructionMethod {
 METHOD_A,
 METHOD_B,
 METHOD_C
};

class NaiveReconstructor {
public:
 std::string reconstructDNA(const DNAInstance &instance,
                            NaiveReconstructionMethod method);

private:
 std::string reconstructDNA_A(const DNAInstance &instance);
 std::string reconstructDNA_B(const DNAInstance &instance);
 std::string reconstructDNA_C(const DNAInstance &instance);
};

int levenshteinDistance(const std::string &s1, const std::string &s2);

#endif //NAIVE_RECONSTRUCTION_H
