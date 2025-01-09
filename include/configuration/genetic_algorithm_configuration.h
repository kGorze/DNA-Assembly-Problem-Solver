//
// Created by konrad_guest on 08/01/2025.
//

#ifndef GENETIC_ALGORITHM_CONFIGURATION_H
#define GENETIC_ALGORITHM_CONFIGURATION_H

#include <memory>
#include "metaheuristics/representation.h"
#include "metaheuristics/selection.h"
#include "metaheuristics/crossover.h"
#include "metaheuristics/mutation.h"
#include "metaheuristics/replacement.h"
#include "metaheuristics/fitness.h"
#include "metaheuristics/stopping_criteria.h"

class GAConfig {
public:
    // Singleton access
    static GAConfig& getInstance();

    // Delete copy constructor and assignment operator
    GAConfig(const GAConfig&) = delete;
    GAConfig& operator=(const GAConfig&) = delete;

    // Configuration getters
    int getPopulationSize() const { return populationSize; }
    int getTournamentSize() const { return tournamentSize; }
    double getMutationRate() const { return mutationRate; }
    int getMaxGenerations() const { return maxGenerations; }
    
    // Component getters
    std::shared_ptr<IRepresentation> getRepresentation() const;
    std::shared_ptr<ISelection> getSelection() const;
    std::shared_ptr<ICrossover> getCrossover(const std::string& type = "order") const;
    std::shared_ptr<IMutation> getMutation() const;
    std::shared_ptr<IReplacement> getReplacement() const;
    std::shared_ptr<IFitness> getFitness() const;
    std::shared_ptr<IStopping> getStopping() const;

    // Custom setters for runtime configuration
    void setPopulationSize(int size) { populationSize = size; }
    void setTournamentSize(int size) { tournamentSize = size; }
    void setMutationRate(double rate) { mutationRate = rate; }
    void setMaxGenerations(int gens) { maxGenerations = gens; }
    void setReplacementRatio(double ratio);

    double getReplacementRatio() const { return m_replacementRatio; }
    std::shared_ptr<IReplacement> getReplacement() {
        return std::make_shared<PartialReplacement>(m_replacementRatio);
    }

private:
    // Private constructor for singleton
    GAConfig();

    // GA parameters
    int populationSize;
    int tournamentSize;
    int maxGenerations;
    double mutationRate = 0.1;
    double m_replacementRatio = 0.7; // default value

};

#endif //GENETIC_ALGORITHM_CONFIGURATION_H
