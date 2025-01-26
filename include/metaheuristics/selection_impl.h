#pragma once

#include "../interfaces/i_selection.h"
#include "../interfaces/i_population_cache.h"
#include "../configuration/genetic_algorithm_configuration.h"
#include "../dna/dna_instance.h"
#include "../utils/logging.h"
#include "utils/random.h"
#include "individual.h"
#include <vector>
#include <memory>
#include <random>
#include <mutex>
#include <algorithm>

class TournamentSelection : public ISelection {
public:
    TournamentSelection(const GAConfig& config, std::shared_ptr<IPopulationCache> cache = nullptr);

    std::vector<std::shared_ptr<Individual>> select(
        const std::vector<std::shared_ptr<Individual>>& population,
        const DNAInstance& instance,
        std::shared_ptr<IFitness> fitness,
        std::shared_ptr<IRepresentation> representation) override;

private:
    const GAConfig& m_config;
    std::shared_ptr<IPopulationCache> m_cache;
}; 