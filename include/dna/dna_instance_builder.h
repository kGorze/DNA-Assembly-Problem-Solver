#ifndef DNA_INSTANCE_BUILDER_H
#define DNA_INSTANCE_BUILDER_H

#include "dna_instance.h"
#include "generator/dna_generator.h"
#include <string>
#include <vector>

class DNAInstanceBuilder {
public:
    DNAInstanceBuilder() = default;
    
    // Builder methods
    DNAInstanceBuilder& setN(int value);
    DNAInstanceBuilder& setK(int value);
    DNAInstanceBuilder& setDeltaK(int value);
    DNAInstanceBuilder& setLNeg(int value);
    DNAInstanceBuilder& setLPoz(int value);
    DNAInstanceBuilder& setRepAllowed(bool value);
    DNAInstanceBuilder& setProbablePositive(int value);
    
    // Build methods
    DNAInstanceBuilder& buildDNA();
    DNAInstanceBuilder& buildSpectrum();
    DNAInstanceBuilder& applyError(IErrorIntroductionStrategy* strategy);
    
    // Get the built instance
    DNAInstance getInstance() const { return instance; }
    
private:
    DNAInstance instance;
    int n = 0;  // DNA length
    ErrorContext errorCtx;
};

#endif // DNA_INSTANCE_BUILDER_H 