#pragma once

#include <string>
#include <vector>
#include "../dna/dna_instance.h"

class SpectrumGenerator {
public:
    /**
     * Generates a spectrum from a DNA sequence.
     * @param dna The DNA sequence to generate spectrum from
     * @param k The base length of oligonucleotides
     * @param deltaK The maximum deviation from k for oligonucleotide length
     * @return Vector of oligonucleotides forming the spectrum
     * @throws std::invalid_argument if parameters are invalid
     */
    static std::vector<std::string> generateSpectrum(const std::string& dna, int k, int deltaK);
    
    /**
     * Generates a spectrum from a DNA instance.
     * @param instance The DNA instance to generate spectrum from
     * @return Vector of oligonucleotides forming the spectrum
     */
    static std::vector<std::string> generateSpectrum(const DNAInstance& instance);
}; 