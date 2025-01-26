#include "../../include/metaheuristics/representation_impl.h"
#include "../../include/utils/logging.h"
#include <set>
#include <sstream>

bool PermutationRepresentation::isValid(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    if (!solution || solution->getGenes().empty()) {
        return false;
    }

    const auto& genes = solution->getGenes();
    const auto& spectrum = instance.getSpectrum();

    // Check if all indices are within bounds
    for (int gene : genes) {
        if (gene < 0 || gene >= static_cast<int>(spectrum.size())) {
            return false;
        }
    }

    // Check for duplicates using a set
    std::set<int> uniqueGenes(genes.begin(), genes.end());
    return uniqueGenes.size() == genes.size();
}

std::string PermutationRepresentation::toString(
    const std::shared_ptr<Individual>& solution,
    [[maybe_unused]] const DNAInstance& instance) const {
    if (!solution) {
        return "";
    }

    std::stringstream ss;
    const auto& genes = solution->getGenes();
    for (size_t i = 0; i < genes.size(); ++i) {
        if (i > 0) {
            ss << " ";
        }
        ss << genes[i];
    }
    return ss.str();
}

std::string PermutationRepresentation::toDNA(
    const std::shared_ptr<Individual>& solution,
    const DNAInstance& instance) const {
    if (!solution) {
        return "";
    }

    std::stringstream ss;
    const auto& genes = solution->getGenes();
    const auto& spectrum = instance.getSpectrum();

    for (int geneIndex : genes) {
        if (geneIndex >= 0 && geneIndex < static_cast<int>(spectrum.size())) {
            ss << spectrum[geneIndex];
        }
    }

    std::string oligo = ss.str();
    if (oligo.size() > static_cast<size_t>(instance.getK())) {
        oligo = oligo.substr(0, instance.getK());
    }
    return oligo;
} 