//
// Created by konrad_guest on 07/01/2025.
//



// ================== SimpleFitness ==================
#include <unordered_set>
#include <iostream>

#include "metaheuristics/fitness.h"


double SimpleFitness::evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // decode
    std::string dna = representation->decodeToDNA(individual, instance);
    int k = instance.getK();
    if(k <= 0) return 0.0;

    // build set from dna's k-mers
    std::unordered_set<std::string> kmerSet;
    for(int i=0; i+(k) <= (int)dna.size(); i++){
        kmerSet.insert(dna.substr(i, k));
    }

    // compare with original instance's spectrum
    const auto &spec = instance.getSpectrum();
    std::unordered_set<std::string> spectrumSet(spec.begin(), spec.end());

    // count how many from kmerSet are in spectrum
    int matches = 0;
    for(auto &ss : kmerSet) {
        if(spectrumSet.count(ss)) {
            matches++;
        }
    }

    return (double)matches;
}

//TODO
// - W kontekście SBH (sekwencjonowania) można zamiast sumy genów liczyć np. liczbę dopasowanych k-merów, minimalizować odległość Levenshteina itp.