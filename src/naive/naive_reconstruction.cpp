//
// Created by konrad_guest on 29/12/2024.
// SMART
#include "naive/naive_reconstruction.h"
#include <algorithm>
#include <random>
#include <sstream>
#include <chrono>

std::string NaiveReconstructor::reconstructDNA(const DNAInstance &instance,
                                               NaiveReconstructionMethod method)
{
    switch(method) {
        case NaiveReconstructionMethod::METHOD_A:
            return reconstructDNA_A(instance);
        case NaiveReconstructionMethod::METHOD_B:
            return reconstructDNA_B(instance);
        case NaiveReconstructionMethod::METHOD_C:
            return reconstructDNA_C(instance);
        default:
            return reconstructDNA_A(instance);
    }
}

std::string NaiveReconstructor::reconstructDNA_A(const DNAInstance &instance)
{
    std::vector<std::string> spectrum = instance.getSpectrum();
    int k = instance.getK();
    if(spectrum.empty()) {
        return "";
    }

    std::sort(spectrum.begin(), spectrum.end());

    std::ostringstream oss;
    oss << spectrum[0];
    for (size_t i = 1; i < spectrum.size(); ++i) {
        oss << spectrum[i].substr(k-1);
    }

    return oss.str();
}

std::string NaiveReconstructor::reconstructDNA_B(const DNAInstance &instance)
{
    std::vector<std::string> spectrum = instance.getSpectrum();
    int k = instance.getK();
    if(spectrum.empty()) {
        return "";
    }

    std::string result = spectrum[0];
    
    for (size_t i = 1; i < spectrum.size(); ++i) {
        const std::string &current = spectrum[i];
        int max_overlap = 0;
        for(int overlap = k-1; overlap > 0; --overlap) {
            if((int)result.size() >= overlap &&
               (int)current.size() >= overlap &&
               result.compare(result.size() - overlap, overlap, current, 0, overlap) == 0) {
                max_overlap = overlap;
                break;
            }
        }
        result += current.substr(max_overlap);
    }

    return result;
}

std::string NaiveReconstructor::reconstructDNA_C(const DNAInstance &instance)
{
    std::vector<std::string> spectrum = instance.getSpectrum();
    int k = instance.getK();
    if(spectrum.empty()) {
        return "";
    }

    auto &rng = RandomGenerator::getInstance().get();
    std::shuffle(spectrum.begin(), spectrum.end(), rng);

    std::ostringstream oss;
    oss << spectrum[0];
    for (size_t i = 1; i < spectrum.size(); ++i) {
        oss << spectrum[i].substr(k-1);
    }
    return oss.str();
}
