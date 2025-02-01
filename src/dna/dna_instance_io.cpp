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

        if (!(file >> n >> k >> deltaK >> lNeg >> lPoz)) {
            throw std::runtime_error("Error reading basic parameters from file");
        }
        
        file >> line; 
        if (line != "0" && line != "1") {
            throw std::runtime_error("Invalid repAllowed value in file");
        }
        repAllowed = (line == "1");
        
        if (!(file >> probablePositive >> size)) {
            throw std::runtime_error("Error reading probability and size from file");
        }
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        try {
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
            if (!std::getline(file, line)) throw std::runtime_error("Error reading DNA sequence");
            instance.setDNA(line);
            
            if (!std::getline(file, line)) throw std::runtime_error("Error reading original DNA sequence");
            instance.setOriginalDNA(line);
            
            if (!std::getline(file, line)) throw std::runtime_error("Error reading target sequence");
            instance.setTargetSequence(line);

            // Read spectrum
            int spectrumSize;
            if (!(file >> spectrumSize)) {
                throw std::runtime_error("Error reading spectrum size");
            }
            file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

            std::vector<std::string> spectrum;
            try {
                spectrum.reserve(spectrumSize);
            } catch (const std::length_error&) {
                throw std::runtime_error("Invalid spectrum size in file");
            }
            
            for (int i = 0; i < spectrumSize; ++i) {
                if (!std::getline(file, line)) {
                    throw std::runtime_error("Error reading spectrum entry");
                }
                spectrum.push_back(line);
            }
            instance.setSpectrum(spectrum);
        } catch (const std::invalid_argument& e) {
            throw std::runtime_error(std::string("Invalid data in file: ") + e.what());
        }

        if (!file.good()) {
            throw std::runtime_error("Error reading from file");
        }

        return instance;
    } catch (const std::runtime_error&) {
        throw;  // Re-throw runtime_error as is
    } catch (const std::exception& e) {
        LOG_ERROR("Error loading instance: " + std::string(e.what()));
        throw std::runtime_error(std::string("Error loading instance: ") + e.what());
    }
} 