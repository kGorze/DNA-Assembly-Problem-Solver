#pragma once

#include "dna/dna_instance.h"
#include "generator/dna_generator.h"
#include "error_introduction.h"
#include <string>
#include <vector>
#include <mutex>
#include <memory>
#include <stdexcept>

class DNAInstanceBuilder {
public:
    DNAInstanceBuilder() = default;
    
    // Builder methods with validation
    DNAInstanceBuilder& setN(int n);
    DNAInstanceBuilder& setK(int k);
    DNAInstanceBuilder& setDeltaK(int deltaK);
    DNAInstanceBuilder& setLNeg(int lNeg);
    DNAInstanceBuilder& setLPoz(int lPoz);
    DNAInstanceBuilder& setRepAllowed(bool repAllowed);
    DNAInstanceBuilder& setProbablePositive(double probablePositive);
    DNAInstanceBuilder& setStartIndex(int startIndex);
    DNAInstanceBuilder& setDNA(const std::string& dna);
    DNAInstanceBuilder& setSpectrum(const std::vector<std::string>& spectrum);
    
    DNAInstanceBuilder& buildDNA();
    DNAInstanceBuilder& buildSpectrum();
    DNAInstanceBuilder& applyError(IErrorIntroductionStrategy* strategy);
    
    // Build method that returns a new instance
    DNAInstance build() const;
    
    // Get a reference to the instance being built
    const DNAInstance& getInstance() const { return m_instance; }
    
private:
    mutable std::mutex m_mutex;
    DNAInstance m_instance;
    DNAGenerator m_generator;
    int m_n = 0;
    int m_k = 0;
    int m_deltaK = 0;
    int m_lPoz = 0;
    int m_lNeg = 0;
    bool m_repAllowed = false;
    double m_probablePositive = 0.0;
    int m_startIndex = 0;
    std::string m_dna;
    std::vector<std::string> m_spectrum;
}; 