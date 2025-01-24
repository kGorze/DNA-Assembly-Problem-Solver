//
// Created by konrad_guest on 07/01/2025.
// SMART

#ifndef SELECTION_H
#define SELECTION_H
#include <vector>
#include <memory>

// Forward declaration
class GAConfig;
class DNAInstance;
class IFitness;
class IRepresentation;
class IPopulationCache;

class ISelection {
public:
    virtual ~ISelection() = default;

    virtual std::vector<std::shared_ptr<std::vector<int>>>
    select(const std::vector<std::shared_ptr<std::vector<int>>> &population,
           const DNAInstance &instance,
           std::shared_ptr<IFitness> fitness,
           std::shared_ptr<IRepresentation> representation) = 0;
};

class TournamentSelection : public ISelection {
public:
    TournamentSelection(GAConfig& config, std::shared_ptr<IPopulationCache> cache);
    std::vector<std::shared_ptr<std::vector<int>>>
    select(const std::vector<std::shared_ptr<std::vector<int>>> &population,
           const DNAInstance &instance,
           std::shared_ptr<IFitness> fitness,
           std::shared_ptr<IRepresentation> representation) override;
    int getTournamentSize() const { return m_tournamentSize; }
private:
    GAConfig& m_config;
    int m_tournamentSize;
    std::shared_ptr<IPopulationCache> m_cache;
};

#endif //SELECTION_H
