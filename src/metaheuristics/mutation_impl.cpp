#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/individual.h"
#include "../../include/utils/logging.h"
#include <memory>
#include <random>
#include <stdexcept>

PointMutation::PointMutation(double mutationRate) {
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        throw std::invalid_argument("Mutation rate must be between 0 and 1");
    }
    m_mutationRate = mutationRate;
    LOG_INFO("PointMutation initialized with rate: " + std::to_string(mutationRate));
}

void PointMutation::mutate(
    std::shared_ptr<Individual>& individual,
    const DNAInstance& instance,
    std::shared_ptr<IRepresentation> representation) {
    
    if (!individual) {
        LOG_ERROR("Null individual provided to mutation operator");
        return;
    }

    if (!representation) {
        LOG_ERROR("Null representation provided to mutation operator");
        return;
    }

    try {
        auto genes = individual->getGenes();
        if (genes.size() < 2) {
            LOG_WARNING("Individual has less than 2 genes, skipping mutation");
            return;
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        std::uniform_int_distribution<size_t> posDis(0, genes.size() - 1);

        bool mutated = false;
        std::vector<int> mutatedGenes = genes;

        // For each position, decide whether to mutate based on mutation rate
        for (size_t i = 0; i < genes.size(); ++i) {
            if (dis(gen) < m_mutationRate) {
                size_t swapPos = posDis(gen);
                while (swapPos == i) {
                    swapPos = posDis(gen);
                }
                
                std::swap(mutatedGenes[i], mutatedGenes[swapPos]);
                mutated = true;
            }
        }

        if (mutated) {
            auto mutatedIndividual = std::make_shared<Individual>(mutatedGenes);
            if (representation->isValid(mutatedIndividual, instance)) {
                individual = mutatedIndividual;
                LOG_DEBUG("Mutation applied successfully");
            } else {
                LOG_WARNING("Mutation produced invalid individual, keeping original");
            }
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Error during mutation: " + std::string(e.what()));
    }
} 