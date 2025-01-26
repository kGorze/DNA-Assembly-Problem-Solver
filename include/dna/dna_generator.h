#pragma once

#include <random>
#include <string>
#include <memory>

class DNAGenerator {
public:
    explicit DNAGenerator(std::shared_ptr<std::mt19937> random);

    void setParameters(int n, int k, int l);

    std::string generateDNA(int length, bool introduceErrors = false);

private:
    std::shared_ptr<std::mt19937> m_random;
    int m_n; // Length of original DNA sequence
    int m_k; // Length of k-mers
    int m_l; // Number of errors
}; 