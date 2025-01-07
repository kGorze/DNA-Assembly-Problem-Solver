//
// Created by konrad_guest on 07/01/2025.
//

#ifndef SELECTION_H
#define SELECTION_H
#include <vector>
#include "metaheuristics/fitness.h"
#include "metaheuristics/representation.h"

/**
 * Interfejs selekcji – przyjmuje populację, zwraca osobniki-rodziców.
 * Zakładamy, że reprezentacja osobnika to wektor double.
 */
class IFitness; // forward declaration

class ISelection {
public:
    virtual ~ISelection() = default;

    virtual std::vector<void*>
    select(const std::vector<void*>                         &population,
                            const DNAInstance               &instance,
                            std::shared_ptr<IFitness>        fitness,
                            std::shared_ptr<IRepresentation> representation) = 0;
};

// ================== Różne rodzaje selekcji ==================

class TournamentSelection : public ISelection {
public:
    TournamentSelection(int tournamentSize);
    std::vector<void*>
    select(const std::vector<void*>                         &population,
                            const DNAInstance               &instance,
                            std::shared_ptr<IFitness>        fitness,
                            std::shared_ptr<IRepresentation> representation) override;
private:
    int m_tournamentSize;
};



// class RouletteSelection : public ISelection {
// public:
//     std::vector<std::vector<double>> select(const std::vector<std::vector<double>> &population,
//                                             const IFitness &fitness) override;
// };
//
// class RankingSelection : public ISelection {
// public:
//     std::vector<std::vector<double>> select(const std::vector<std::vector<double>> &population,
//                                             const IFitness &fitness) override;
// };
//
// class ElitistSelection : public ISelection {
// public:
//     ElitistSelection(int eliteCount);
//     std::vector<std::vector<double>> select(const std::vector<std::vector<double>> &population,
//                                             const IFitness &fitness) override;
// private:
//     int m_eliteCount;
// };

#endif //SELECTION_H
