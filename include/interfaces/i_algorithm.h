#pragma once

#include "dna/dna_instance.h"
#include <string>

class IAlgorithm {
public:
    virtual ~IAlgorithm() = default;
    virtual std::string run(const DNAInstance& instance) = 0;
}; 