//
// Created by konrad_guest on 07/01/2025.
// SMART

#include "../../include/interfaces/i_representation.h"
#include "../../include/metaheuristics/representation_impl.h"
#include "../../include/utils/logging.h"
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

/************************************************
 *  DirectDNARepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
DirectDNARepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    int n = instance.getN();
    if (n <= 0) n = fallbackN; 

    std::mt19937 rng(std::random_device{}());
    static const char letters[] = "ACGT";

    for (int i = 0; i < popSize; i++) {
        auto individual = std::make_shared<std::vector<int>>();
        individual->reserve(n);

        for (int j = 0; j < n; j++) {
            int base = rng() % 4; // 0..3
            individual->push_back(base);
        }
        population.push_back(individual);
    }

    return population;
}

bool DirectDNARepresentation::isValid(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    if (!solution || solution->empty()) return false;
    for (int val : *solution) {
        if (val < 0 || val > 3) return false;
    }
    return true;
}

std::string DirectDNARepresentation::toString(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    static const char letters[] = {'A','C','G','T'};
    if (!solution) return "";

    std::string result;
    result.reserve(solution->size());
    for (auto val : *solution) {
        int idx = std::max(0, std::min(val, 3));
        result.push_back(letters[idx]);
    }
    return result;
}

std::vector<char> DirectDNARepresentation::toDNA(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    std::string str = toString(solution, instance);
    return std::vector<char>(str.begin(), str.end());
}

/************************************************
 *  PermutationRepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
PermutationRepresentation::initializePopulation(
    int populationSize,
    const DNAInstance& instance) 
{
    const size_t spectrumSize = instance.getSpectrum().size();
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(populationSize);

    for (int i = 0; i < populationSize; ++i) {
        auto individual = std::make_shared<std::vector<int>>(spectrumSize);
        std::iota(individual->begin(), individual->end(), 0);
        std::random_shuffle(individual->begin(), individual->end());
        
        if (isValid(individual, instance)) {
            population.push_back(individual);
        } else {
            --i;
        }
    }
    return population;
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    if (!solution || solution->empty()) return false;
    const size_t spectrumSize = instance.getSpectrum().size();
    std::vector<bool> used(spectrumSize, false);
    
    for (int val : *solution) {
        if (val < 0 || val >= static_cast<int>(spectrumSize) || used[val]) return false;
        used[val] = true;
    }
    return true;
}

std::string PermutationRepresentation::toString(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    if (!solution || solution->empty()) {
        LOG_WARNING("Solution is nullptr or empty");
        return "";
    }
    
    const auto& spec = instance.getSpectrum();
    int k = instance.getK();
    int n = instance.getN();
    
    if (k <= 0 || n <= 0) {
        LOG_WARNING("Invalid k or n values");
        return "";
    }

    int firstIndex = (*solution)[0];
    if (firstIndex < 0 || firstIndex >= static_cast<int>(spec.size())) {
        LOG_WARNING("Invalid first index in solution");
        return "";
    }

    std::string result = spec[firstIndex];

    for (size_t i = 1; i < solution->size(); i++) {
        int fragIdx = (*solution)[i];
        if (fragIdx < 0 || fragIdx >= static_cast<int>(spec.size())) continue;

        const std::string& frag = spec[fragIdx];
        if (static_cast<int>(frag.size()) < k - 2) continue;
        result.push_back(frag.back());
    }

    if (static_cast<int>(result.size()) > n) {
        result.resize(n);
    } else while (static_cast<int>(result.size()) < n) {
        result.push_back('N');
    }

    return result;
}

std::vector<char> PermutationRepresentation::toDNA(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    std::string str = toString(solution, instance);
    return std::vector<char>(str.begin(), str.end());
}

/************************************************
 *  GraphPathRepresentation
 ************************************************/
std::vector<std::shared_ptr<std::vector<int>>>
GraphPathRepresentation::initializePopulation(int popSize, const DNAInstance &instance)
{
    std::vector<std::shared_ptr<std::vector<int>>> population;
    population.reserve(popSize);

    int k = instance.getK();
    int n = instance.getN();
    int requiredLength = n - k + 1;
    
    for(int i = 0; i < popSize; i++) {
        auto path = std::make_shared<std::vector<int>>();
        path->reserve(requiredLength);
        
        for(int idx = 0; idx < static_cast<int>(instance.getSpectrum().size()); idx++) {
            path->push_back(idx);
        }
        
        std::shuffle(path->begin(), path->end(), std::mt19937(std::random_device{}()));
        path->resize(requiredLength);
        
        population.push_back(path);
    }

    return population;
}

bool GraphPathRepresentation::isValid(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    if (!solution || solution->empty()) return false;
    
    int k = instance.getK();
    int n = instance.getN();
    int requiredLength = n - k + 1;
    
    if (static_cast<int>(solution->size()) != requiredLength) return false;
    
    const auto& spec = instance.getSpectrum();
    for (int val : *solution) {
        if (val < 0 || val >= static_cast<int>(spec.size())) return false;
    }
    
    return true;
}

std::string GraphPathRepresentation::toString(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    if (!solution || solution->empty()) return "";
    
    const auto& spec = instance.getSpectrum();
    if (spec.empty()) return "";

    int k = instance.getK();
    int n = instance.getN();
    
    int firstIndex = (*solution)[0];
    if (firstIndex < 0 || firstIndex >= static_cast<int>(spec.size())) return "";

    std::string result = spec[firstIndex];

    for (size_t i = 1; i < solution->size(); i++) {
        int fragIdx = (*solution)[i];
        if (fragIdx < 0 || fragIdx >= static_cast<int>(spec.size())) continue;

        const std::string& frag = spec[fragIdx];
        if (static_cast<int>(frag.size()) < k - 2) continue;
        result.push_back(frag.back());
    }

    if (static_cast<int>(result.size()) != n) {
        LOG_WARNING("Invalid DNA length");
        return "";
    }

    return result;
}

std::vector<char> GraphPathRepresentation::toDNA(
    const std::shared_ptr<std::vector<int>>& solution,
    const DNAInstance& instance) const 
{
    std::string str = toString(solution, instance);
    return std::vector<char>(str.begin(), str.end());
}
