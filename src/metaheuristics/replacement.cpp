//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "metaheuristics/replacement.h"
#include "../include/metaheuristics/replacement_impl.h"
#include "../include/utils/logging.h"
#include <algorithm>
#include <iostream>

std::vector<std::shared_ptr<std::vector<int>>>
FullReplacement::replace(const std::vector<std::shared_ptr<std::vector<int>>> &oldPop,
                         const std::vector<std::shared_ptr<std::vector<int>>> &offspring,
                         const DNAInstance &instance,
                         std::shared_ptr<IFitness> fitness,
                         std::shared_ptr<IRepresentation> representation)
{
    // trivial: discard oldPop, return offspring
    return offspring;
}

// =========== PartialReplacement ===========

PartialReplacement::PartialReplacement(double replacementRatio, std::shared_ptr<IPopulationCache> cache) 
    : m_replacementRatio(replacementRatio)
    , m_fitnessCache(cache) {}

PartialReplacement::~PartialReplacement() = default;

std::vector<std::shared_ptr<std::vector<int>>> 
PartialReplacement::replace(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
    const std::vector<double>& populationFitness,
    const std::vector<double>& offspringFitness,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (population.empty() || offspring.empty()) {
        return population;
    }

    // Calculate how many offspring to include
    size_t numToReplace = static_cast<size_t>(population.size() * m_replacementRatio);
    numToReplace = std::min(numToReplace, offspring.size());

    // Create combined population
    std::vector<std::pair<double, std::shared_ptr<std::vector<int>>>> combined;
    combined.reserve(population.size() + offspring.size());

    // Add population with their fitness
    for (size_t i = 0; i < population.size(); ++i) {
        combined.emplace_back(populationFitness[i], population[i]);
    }

    // Add offspring with their fitness
    for (size_t i = 0; i < offspring.size(); ++i) {
        combined.emplace_back(offspringFitness[i], offspring[i]);
    }

    // Sort by fitness in descending order
    std::sort(combined.begin(), combined.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // Take the best individuals
    std::vector<std::shared_ptr<std::vector<int>>> newPopulation;
    newPopulation.reserve(population.size());

    for (size_t i = 0; i < population.size(); ++i) {
        newPopulation.push_back(combined[i].second);
    }

    return newPopulation;
}

std::vector<std::shared_ptr<std::vector<int>>> 
ElitistReplacement::replace(
    const std::vector<std::shared_ptr<std::vector<int>>>& population,
    const std::vector<std::shared_ptr<std::vector<int>>>& offspring,
    const std::vector<double>& populationFitness,
    const std::vector<double>& offspringFitness,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation)
{
    if (population.empty() || offspring.empty()) {
        return population;
    }

    // Find the best individual from the current population
    auto maxPopIt = std::max_element(populationFitness.begin(), populationFitness.end());
    size_t bestPopIdx = std::distance(populationFitness.begin(), maxPopIt);
    auto bestFromPop = population[bestPopIdx];
    double bestPopFitness = *maxPopIt;

    // Create new population starting with the best from the current population
    std::vector<std::shared_ptr<std::vector<int>>> newPopulation;
    newPopulation.reserve(population.size());
    newPopulation.push_back(bestFromPop);

    // Add the best offspring until the population is full
    std::vector<std::pair<double, std::shared_ptr<std::vector<int>>>> sortedOffspring;
    sortedOffspring.reserve(offspring.size());
    for (size_t i = 0; i < offspring.size(); ++i) {
        sortedOffspring.emplace_back(offspringFitness[i], offspring[i]);
    }

    std::sort(sortedOffspring.begin(), sortedOffspring.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    for (size_t i = 0; i < population.size() - 1 && i < sortedOffspring.size(); ++i) {
        newPopulation.push_back(sortedOffspring[i].second);
    }

    // If we still need more individuals, take them from the original population
    while (newPopulation.size() < population.size()) {
        newPopulation.push_back(population[newPopulation.size() - 1]);
    }

    return newPopulation;
}
