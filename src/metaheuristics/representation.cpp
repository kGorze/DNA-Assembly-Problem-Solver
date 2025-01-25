//
// Created by konrad_guest on 07/01/2025.
// SMART

#include "../../include/interfaces/i_representation.h"
#include "../../include/metaheuristics/representation_impl.h"
#include "../../include/utils/logging.h"
#include "../../include/metaheuristics/individual.h"
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
std::vector<std::shared_ptr<Individual>>
PermutationRepresentation::initializePopulation(int popSize, const DNAInstance& instance)
{
    if (popSize <= 0) {
        LOG_ERROR("Invalid population size: " + std::to_string(popSize));
        return {};
    }

    std::vector<std::shared_ptr<Individual>> population;
    population.reserve(popSize);

    int validCount = 0;
    int invalidCount = 0;

    for (int i = 0; i < popSize; i++) {
        try {
            auto individual = std::make_shared<Individual>();
            if (initializeIndividual(*individual, instance)) {
                population.push_back(individual);
                validCount++;
            } else {
                LOG_WARNING("Failed to initialize individual " + std::to_string(i));
                invalidCount++;
            }
        } catch (const std::exception& e) {
            LOG_ERROR("Exception while initializing individual: " + std::string(e.what()));
            invalidCount++;
        }
    }

    LOG_INFO("Population initialization complete:");
    LOG_INFO("  Valid individuals: " + std::to_string(validCount));
    LOG_INFO("  Invalid individuals: " + std::to_string(invalidCount));

    return population;
}

bool PermutationRepresentation::initializeIndividual(Individual& individual, const DNAInstance& instance)
{
    try {
        auto genes = generateRandomSolution(instance);
        if (!validateSolution(genes, instance)) {
            LOG_ERROR("Generated invalid solution");
            return false;
        }
        individual.setGenes(std::move(genes));
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to initialize individual: " + std::string(e.what()));
        return false;
    }
}

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to isValid");
        return false;
    }

    try {
        const auto& genes = solution->getGenes();
        return validateSolution(genes, instance);
    } catch (const std::exception& e) {
        LOG_ERROR("Exception in isValid: " + std::string(e.what()));
        return false;
    }
}

std::string PermutationRepresentation::toString(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to toString");
        return "";
    }

    try {
        return solution->toString();
    } catch (const std::exception& e) {
        LOG_ERROR("Exception in toString: " + std::string(e.what()));
        return "";
    }
}

std::vector<char> PermutationRepresentation::toDNA(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const
{
    if (!solution) {
        LOG_ERROR("Null solution provided to toDNA");
        return {};
    }

    try {
        std::string str = toString(solution, instance);
        return std::vector<char>(str.begin(), str.end());
    } catch (const std::exception& e) {
        LOG_ERROR("Exception in toDNA: " + std::string(e.what()));
        return {};
    }
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
