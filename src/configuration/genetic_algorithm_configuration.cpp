//
// Created by konrad_guest on 08/01/2025.
//
#include "configuration/genetic_algorithm_configuration.h"

GAConfig::GAConfig() {
    // Default parameters
    populationSize = 100;
    tournamentSize = 3;
    mutationRate = 0.1;
    maxGenerations = 200;
    
}

GAConfig& GAConfig::getInstance() {
    static GAConfig instance;
    return instance;
}

std::shared_ptr<IRepresentation> GAConfig::getRepresentation() const {
    return std::make_shared<PermutationRepresentation>();
}

std::shared_ptr<ISelection> GAConfig::getSelection() const {
    return std::make_shared<TournamentSelection>(tournamentSize);
}

std::shared_ptr<ICrossover> GAConfig::getCrossover(const std::string& type) const {
    if (type == "order") {
        return std::make_shared<OrderCrossover>();
    } else if (type == "edge") {
        return std::make_shared<EdgeRecombination>();
    } else if (type == "pmx") {
        return std::make_shared<PMXCrossover>();
    } else if (type == "onepoint") {
        return std::make_shared<OnePointCrossover>();
    }
    // Default to OrderCrossover if type is not recognized
    return std::make_shared<OrderCrossover>();
}

std::shared_ptr<IMutation> GAConfig::getMutation() const {
    auto mutation = std::make_shared<PointMutation>(mutationRate);
    return mutation;
}

std::shared_ptr<IReplacement> GAConfig::getReplacement() const {
    return std::make_shared<PartialReplacement>();
}

std::shared_ptr<IFitness> GAConfig::getFitness() const {
    return std::make_shared<BetterFitness>();
}

std::shared_ptr<IStopping> GAConfig::getStopping() const {
    return std::make_shared<MaxGenerationsStopping>(maxGenerations);
}

void GAConfig::setReplacementRatio(double ratio) {
    if (ratio < 0.0) ratio = 0.0;
    if (ratio > 1.0) ratio = 1.0;
    m_replacementRatio = ratio;
}