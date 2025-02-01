#include "dna/dna_instance_io.h"
#include "utils/logging.h"
#include <fstream>
#include <stdexcept>

bool InstanceIO::saveInstance(const DNAInstance& instance, const std::string& filename) {
    try {
        std::ofstream file(filename);
        if (!file) {
            LOG_ERROR("Could not open file for writing: " + filename);
            throw std::runtime_error("Could not open file for writing: " + filename);
        }

        // Write basic parameters with high precision for floating point values
        file.precision(17);  // Use maximum precision for double values
        file << instance.getN() << "\n";
        file << instance.getK() << "\n";
        file << instance.getDeltaK() << "\n";
        file << instance.getLNeg() << "\n";
        file << instance.getLPoz() << "\n";
        file << (instance.isRepAllowed() ? "1" : "0") << "\n";
        file << instance.getProbablePositive() << "\n";
        file << instance.getSize() << "\n";

        // Write DNA sequences
        file << instance.getDNA() << "\n";
        file << instance.getOriginalDNA() << "\n";
        file << instance.getTargetSequence() << "\n";

        // Write spectrum
        const auto& spectrum = instance.getSpectrum();
        file << spectrum.size() << "\n";
        for (const auto& kmer : spectrum) {
            file << kmer << "\n";
        }

        return file.good();
    } catch (const std::exception& e) {
        LOG_ERROR("Error saving instance: " + std::string(e.what()));
        throw;
    }
}

DNAInstance InstanceIO::loadInstance(const std::string& filename) {
    try {
        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file for reading: " + filename);
        }

        DNAInstance instance;
        std::string line;

        // Read basic parameters
        int n, k, deltaK, lNeg, lPoz, size;
        bool repAllowed;
        double probablePositive;

        file >> n >> k >> deltaK >> lNeg >> lPoz;
        file >> line; repAllowed = (line == "1");
        file >> probablePositive;
        file >> size;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');  // Skip to next line

        // Set basic parameters
        instance.setN(n);
        instance.setK(k);
        instance.setDeltaK(deltaK);
        instance.setLNeg(lNeg);
        instance.setLPoz(lPoz);
        instance.setRepAllowed(repAllowed);
        instance.setProbablePositive(probablePositive);
        instance.setSize(size);

        // Read DNA sequences
        std::getline(file, line); instance.setDNA(line);
        std::getline(file, line); instance.setOriginalDNA(line);
        std::getline(file, line); instance.setTargetSequence(line);

        // Read spectrum
        int spectrumSize;
        file >> spectrumSize;
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        std::vector<std::string> spectrum;
        spectrum.reserve(spectrumSize);
        for (int i = 0; i < spectrumSize; ++i) {
            std::getline(file, line);
            spectrum.push_back(line);
        }
        instance.setSpectrum(spectrum);

        if (!file.good()) {
            throw std::runtime_error("Error reading from file");
        }

        return instance;
    } catch (const std::exception& e) {
        LOG_ERROR("Error loading instance: " + std::string(e.what()));
        throw;
    }
} 