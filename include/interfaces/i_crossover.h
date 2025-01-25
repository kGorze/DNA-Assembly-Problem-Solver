#pragma once

#include <memory>
#include <vector>
#include "../dna/dna_instance.h"
#include "i_representation.h"

class ICrossover {
public:
    virtual ~ICrossover() = default;

    virtual std::vector<std::shared_ptr<std::vector<int>>> crossover(
        const std::vector<std::shared_ptr<std::vector<int>>>& parents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation) = 0;
}; 