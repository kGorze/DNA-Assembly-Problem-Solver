//
// Created by konrad_guest on 07/01/2025.
// SMART
#include "../include/metaheuristics/replacement_impl.h"
#include "../include/utils/logging.h"
#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>

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
        LOG_WARNING("Empty population or offspring in replacement");
        return population;
    }

    // Calculate number of individuals to replace
    size_t numToReplace = static_cast<size_t>(m_replacementRatio * population.size());
    numToReplace = std::min(numToReplace, offspring.size());

    // Combine populations and their fitness values
    std::vector<std::pair<double, std::shared_ptr<std::vector<int>>>> combined;
    combined.reserve(population.size() + offspring.size());

    for (size_t i = 0; i < population.size(); ++i) {
        combined.emplace_back(populationFitness[i], population[i]);
    }
    for (size_t i = 0; i < offspring.size(); ++i) {
        combined.emplace_back(offspringFitness[i], offspring[i]);
    }

    // Sort by fitness in descending order
    std::sort(combined.begin(), combined.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // Take the best individuals
    std::vector<std::shared_ptr<std::vector<int>>> newPopulation;
    newPopulation.reserve(population.size());

    // Keep the best from the original population
    size_t keepFromPopulation = population.size() - numToReplace;
    for (size_t i = 0; i < keepFromPopulation && i < combined.size(); ++i) {
        newPopulation.push_back(combined[i].second);
    }

    // Add the best offspring
    for (size_t i = 0; i < numToReplace && keepFromPopulation + i < combined.size(); ++i) {
        newPopulation.push_back(combined[keepFromPopulation + i].second);
    }

    std::stringstream ss;
    ss << "Replaced " << numToReplace << " individuals";
    LOG_INFO(ss.str());

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
        LOG_WARNING("Empty population or offspring in replacement");
        return population;
    }

    // Find the best individual from current population
    auto bestPopIt = std::max_element(populationFitness.begin(), populationFitness.end());
    size_t bestPopIdx = std::distance(populationFitness.begin(), bestPopIt);

    // Create new population starting with the best from current population
    std::vector<std::shared_ptr<std::vector<int>>> newPopulation;
    newPopulation.reserve(population.size());
    newPopulation.push_back(population[bestPopIdx]);

    // Sort offspring by fitness
    std::vector<std::pair<double, std::shared_ptr<std::vector<int>>>> sortedOffspring;
    sortedOffspring.reserve(offspring.size());
    for (size_t i = 0; i < offspring.size(); ++i) {
        sortedOffspring.emplace_back(offspringFitness[i], offspring[i]);
    }
    std::sort(sortedOffspring.begin(), sortedOffspring.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    // Fill rest of population with best offspring
    for (size_t i = 0; i < population.size() - 1 && i < sortedOffspring.size(); ++i) {
        newPopulation.push_back(sortedOffspring[i].second);
    }

    // If we still need more individuals, take them from original population
    while (newPopulation.size() < population.size()) {
        size_t idx = newPopulation.size() - 1;
        if (idx < population.size()) {
            newPopulation.push_back(population[idx]);
        }
    }

    std::stringstream ss;
    ss << "Elitist replacement complete. Best fitness: " << *bestPopIt;
    LOG_INFO(ss.str());

    return newPopulation;
}
