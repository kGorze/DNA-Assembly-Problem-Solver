#pragma once

#include <string>
#include <vector>
#include <memory>
#include <random>

class SpectrumGenerator {
public:
    explicit SpectrumGenerator(std::shared_ptr<std::mt19937> random);

    std::vector<std::string> generateSpectrum(const std::string& dna, int k, int errorCount = 0);

private:
    void generateErrorKmers(const std::string& kmer, int errorsLeft, std::vector<std::string>& result);
    std::shared_ptr<std::mt19937> m_random;
}; 