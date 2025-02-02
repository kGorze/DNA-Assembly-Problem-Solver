//test already in the suite.
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <random>

class SpectrumGenerator {
private:
    std::shared_ptr<std::mt19937> m_random;
    void generateErrorKmers(const std::string& kmer, int errorsLeft, std::vector<std::string>& result);

public:
    explicit SpectrumGenerator(std::shared_ptr<std::mt19937> random = nullptr);
    std::vector<std::string> generateSpectrum(const std::string& dna, int k, int errorCount);
}; 