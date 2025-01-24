#include "dna/dna_instance_io.h"
#include "utils/logging.h"
#include <fstream>
#include <sstream>

bool InstanceIO::loadInstance(const std::string& filename, DNAInstance& instance) {
    std::ifstream in(filename);
    if (!in) {
        LOG_ERROR("Could not open file for reading: " + filename);
        return false;
    }

    // Read parameters
    int n, k, deltaK, lNeg, lPoz, probablePositive;
    bool repAllowed;
    in >> n >> k >> deltaK >> lNeg >> lPoz >> repAllowed >> probablePositive;

    instance.setN(n);
    instance.setK(k);
    instance.setDeltaK(deltaK);
    instance.setLNeg(lNeg);
    instance.setLPoz(lPoz);
    instance.setRepAllowed(repAllowed);
    instance.setProbablePositive(probablePositive);

    // Read DNA sequence
    std::string dna;
    in >> dna;
    instance.setDNA(dna);

    // Read spectrum
    int spectrumSize;
    in >> spectrumSize;
    std::vector<std::string> spectrum;
    spectrum.reserve(spectrumSize);

    std::string oligo;
    for (int i = 0; i < spectrumSize; ++i) {
        in >> oligo;
        spectrum.push_back(oligo);
    }
    instance.setSpectrum(spectrum);

    return true;
}

bool InstanceIO::saveInstance(const DNAInstance& instance, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        LOG_ERROR("Could not open file for writing: " + filename);
        return false;
    }

    // Write parameters
    out << instance.getN() << " "
        << instance.getK() << " "
        << instance.getDeltaK() << " "
        << instance.getLNeg() << " "
        << instance.getLPoz() << " "
        << instance.isRepAllowed() << " "
        << instance.getProbablePositive() << "\n";

    // Write DNA sequence
    out << instance.getDNA() << "\n";

    // Write spectrum
    const auto& spectrum = instance.getSpectrum();
    out << spectrum.size() << "\n";
    for (const auto& oligo : spectrum) {
        out << oligo << "\n";
    }

    return true;
} 