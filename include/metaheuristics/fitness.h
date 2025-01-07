//
// Created by konrad_guest on 07/01/2025.
//

#ifndef FITNESS_H
#define FITNESS_H

#include "metaheuristics/representation.h"
#include "generator/dna_generator.h"

#include <vector>

class IFitness {
public:
    virtual ~IFitness() = default;
    /**
     * Ocena przystosowania pojedynczego osobnika
     */
    virtual double evaluate(void* individual,
                        const DNAInstance &instance,
                        std::shared_ptr<IRepresentation> representation) const = 0;
};

/**
 * Przykładowa implementacja – do rozbudowy
 */
class SimpleFitness : public IFitness {
public:
    double evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const override;
};

#endif //FITNESS_H
