#include "dna/dna_instance_io.h"
#include "utils/logging.h"
#include <fstream>
#include <stdexcept>

bool InstanceIO::saveInstance(const DNAInstance& instance, const std::string& filename) {
    try {
        std::ofstream file(filename);
        if (!file) {
            LOG_ERROR("Could not open file for writing: " + filename);
            return false;
        }

        // Write basic parameters
        file << instance.getN() << "\n";
        file << instance.getK() << "\n";
        file << instance.getDeltaK() << "\n";
        file << instance.getLNeg() << "\n";
        file << instance.getLPoz() << "\n";
        file << instance.isRepAllowed() << "\n";
        file << instance.getProbablePositive() << "\n";

        // Write DNA sequence and target sequence
        file << instance.getDNA() << "\n";
        file << instance.getTargetSequence() << "\n";

        // Write spectrum
        const auto& spectrum = instance.getSpectrum();
        file << spectrum.size() << "\n";
        for (const auto& kmer : spectrum) {
            file << kmer << "\n";
        }

        // Write start index
        file << instance.getStartIndex() << "\n";

        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Error saving instance: " + std::string(e.what()));
        return false;
    }
}

DNAInstance InstanceIO::loadInstance(const std::string& filename) {
    try {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file for reading: " + filename);
        }

        DNAInstance instance;
        int value;

        // Read basic parameters
        file >> value;
        instance.setN(value);
        file >> value;
        instance.setK(value);
        file >> value;
        instance.setDeltaK(value);
        file >> value;
        instance.setLNeg(value);
        file >> value;
        instance.setLPoz(value);
        
        bool repAllowed;
        file >> repAllowed;
        instance.setRepAllowed(repAllowed);
        
        file >> value;
        instance.setProbablePositive(value);

        // Read DNA sequence and target sequence
        std::string dna, target;
        file >> dna;
        file >> target;
        instance.setDNA(dna);
        instance.setTargetSequence(target);

        // Read spectrum
        int spectrumSize;
        file >> spectrumSize;
        std::vector<std::string> spectrum;
        spectrum.reserve(spectrumSize);
        
        std::string kmer;
        for (int i = 0; i < spectrumSize; ++i) {
            file >> kmer;
            spectrum.push_back(kmer);
        }
        instance.setSpectrum(spectrum);

        // Read start index
        int startIndex;
        file >> startIndex;
        instance.setStartIndex(startIndex);

        return instance;
    } catch (const std::exception& e) {
        LOG_ERROR("Error loading instance: " + std::string(e.what()));
        throw;
    }
} 