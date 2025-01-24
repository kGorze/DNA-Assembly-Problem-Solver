#pragma once

#include "../interfaces/i_selection.h"
#include "../dna/dna_instance.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include <vector>
#include <memory>
#include <random>

class TournamentSelection : public ISelection {
public:
    TournamentSelection(GAConfig& config, std::shared_ptr<IPopulationCache> cache) 
        : m_config(config), m_tournamentSize(config.getTournamentSize()), m_cache(cache) {}
    
    std::vector<std::shared_ptr<std::vector<int>>> select(
        const std::vector<std::shared_ptr<std::vector<int>>>& population,
        const std::vector<double>& fitness,
        size_t numParents,
        const DNAInstance& instance,
        std::shared_ptr<IRepresentation> representation
    ) override;

private:
    GAConfig& m_config;
    int m_tournamentSize;
    std::shared_ptr<IPopulationCache> m_cache;
}; 