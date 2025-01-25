#ifndef DNA_INSTANCE_IO_H
#define DNA_INSTANCE_IO_H

#include "dna_instance.h"
#include <string>
#include <mutex>
#include <stdexcept>

class InstanceIO {
public:
    static bool loadInstance(const std::string& filename, DNAInstance& instance) {
        if (filename.empty()) {
            throw std::invalid_argument("Filename cannot be empty");
        }
        
        std::lock_guard<std::mutex> lock(s_mutex);
        try {
            std::ifstream in(filename);
            if (!in) {
                LOG_ERROR("Could not open file for reading: " + filename);
                return false;
            }

            // Read parameters
            int n, k, deltaK, lNeg, lPoz, probablePositive;
            bool repAllowed;
            in >> n >> k >> deltaK >> lNeg >> lPoz >> repAllowed >> probablePositive;

            if (n <= 0 || k <= 0 || deltaK < 0 || lNeg < 0 || lPoz < 0) {
                LOG_ERROR("Invalid parameters in file: " + filename);
                return false;
            }

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
            if (dna.empty()) {
                LOG_ERROR("Empty DNA sequence in file: " + filename);
                return false;
            }
            instance.setDNA(dna);

            // Read spectrum
            int spectrumSize;
            in >> spectrumSize;
            if (spectrumSize <= 0) {
                LOG_ERROR("Invalid spectrum size in file: " + filename);
                return false;
            }

            std::vector<std::string> spectrum;
            spectrum.reserve(spectrumSize);

            std::string oligo;
            for (int i = 0; i < spectrumSize; ++i) {
                in >> oligo;
                if (oligo.empty()) {
                    LOG_ERROR("Empty oligonucleotide in spectrum at position " + std::to_string(i));
                    return false;
                }
                spectrum.push_back(oligo);
            }
            instance.setSpectrum(spectrum);

            return true;
        } catch (const std::exception& e) {
            LOG_ERROR("Error loading instance from file: " + std::string(e.what()));
            return false;
        }
    }

    static bool saveInstance(const DNAInstance& instance, const std::string& filename) {
        if (filename.empty()) {
            throw std::invalid_argument("Filename cannot be empty");
        }
        
        std::lock_guard<std::mutex> lock(s_mutex);
        try {
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
            const std::string dna = instance.getDNA();
            if (dna.empty()) {
                LOG_ERROR("Cannot save empty DNA sequence");
                return false;
            }
            out << dna << "\n";

            // Write spectrum
            const auto spectrum = instance.getSpectrum();
            if (spectrum.empty()) {
                LOG_ERROR("Cannot save empty spectrum");
                return false;
            }
            out << spectrum.size() << "\n";
            for (const auto& oligo : spectrum) {
                if (oligo.empty()) {
                    LOG_ERROR("Cannot save empty oligonucleotide in spectrum");
                    return false;
                }
                out << oligo << "\n";
            }

            return true;
        } catch (const std::exception& e) {
            LOG_ERROR("Error saving instance to file: " + std::string(e.what()));
            return false;
        }
    }

private:
    static std::mutex s_mutex;  // Mutex for thread-safe file operations
};

std::mutex InstanceIO::s_mutex;  // Definition of static mutex

#endif // DNA_INSTANCE_IO_H 