#pragma once

#include "../interfaces/i_selection.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include <vector>
#include <memory>
#include <random>

class TournamentSelection : public ISelection {
public:
    TournamentSelection(GAConfig& config, std::shared_ptr<IPopulationCache> cache);
    
    std::vector<std::shared_ptr<std::vector<int>>> select(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation
    ) override;

private:
    GAConfig& m_config;
    int m_tournamentSize;
    std::shared_ptr<IPopulationCache> m_cache;
}; 