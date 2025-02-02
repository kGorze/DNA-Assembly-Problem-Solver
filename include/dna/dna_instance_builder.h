//test already in the suite.
#pragma once

#include "dna/dna_instance.h"
#include "dna/error_introduction.h"
#include "generator/dna_generator.h"
#include <string>
#include <vector>
#include <mutex>
#include <memory>
#include <stdexcept>

class DNAInstanceBuilder {
public:
    explicit DNAInstanceBuilder(std::shared_ptr<DNAGenerator> generator);
    
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
    bool validateState(const std::string& context) const;
    mutable std::mutex m_mutex;
    DNAInstance m_instance;
    std::shared_ptr<DNAGenerator> m_generator;
}; 