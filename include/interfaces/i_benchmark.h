#pragma once

#include "../dna/dna_instance.h"
#include <string>

class IBenchmark {
public:
    virtual ~IBenchmark() = default;
    virtual void runBenchmark(const DNAInstance& instance) = 0;
    virtual std::string getName() const = 0;
}; 