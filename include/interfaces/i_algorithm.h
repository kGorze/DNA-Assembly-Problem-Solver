#pragma once

#include "dna/dna_instance.h"

class IAlgorithm {
public:
    virtual ~IAlgorithm() = default;
    virtual void run(const DNAInstance& instance) = 0;
}; 