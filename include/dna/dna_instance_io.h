#ifndef DNA_INSTANCE_IO_H
#define DNA_INSTANCE_IO_H

#include "dna_instance.h"
#include <string>
#include <mutex>
#include <stdexcept>
#include <fstream>
#include <sstream>

class InstanceIO {
public:
    static bool loadInstance(const std::string& filename, DNAInstance& instance) {
        if (filename.empty()) {
            throw std::invalid_argument("Filename cannot be empty");
        }
        
        static std::mutex s_mutex;
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
            
            if (in.fail()) {
                LOG_ERROR("Failed to read parameters from file: " + filename);
                return false;
            }
            
            // Validate parameters
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
            
            if (in.fail() || dna.empty()) {
                LOG_ERROR("Failed to read DNA sequence from file: " + filename);
                return false;
            }
            
            instance.setDNA(dna);
            
            // Read spectrum
            int spectrumSize;
            in >> spectrumSize;
            
            if (in.fail() || spectrumSize < 0) {
                LOG_ERROR("Failed to read spectrum size from file: " + filename);
                return false;
            }
            
            std::vector<std::string> spectrum;
            spectrum.reserve(spectrumSize);
            
            std::string oligo;
            for (int i = 0; i < spectrumSize; ++i) {
                in >> oligo;
                if (in.fail() || oligo.empty()) {
                    LOG_ERROR("Failed to read oligonucleotide from file: " + filename);
                    return false;
                }
                spectrum.push_back(oligo);
            }
            
            instance.setSpectrum(spectrum);
            
            // Find and set start index
            std::string startFrag = dna.substr(0, k);
            int startIdx = -1;
            for (int i = 0; i < static_cast<int>(spectrum.size()); i++) {
                if (spectrum[i] == startFrag) {
                    startIdx = i;
                    break;
                }
            }
            instance.setStartIndex(startIdx);
            
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
        
        static std::mutex s_mutex;
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
            
            if (out.fail()) {
                LOG_ERROR("Failed to write parameters to file: " + filename);
                return false;
            }
            
            // Write DNA sequence
            const std::string dna = instance.getDNA();
            if (dna.empty()) {
                LOG_ERROR("Cannot save empty DNA sequence to file: " + filename);
                return false;
            }
            out << dna << "\n";
            
            // Write spectrum
            const auto& spectrum = instance.getSpectrum();
            if (spectrum.empty()) {
                LOG_ERROR("Cannot save empty spectrum to file: " + filename);
                return false;
            }
            
            out << spectrum.size() << "\n";
            for (const auto& oligo : spectrum) {
                out << oligo << "\n";
            }
            
            if (out.fail()) {
                LOG_ERROR("Failed to write spectrum to file: " + filename);
                return false;
            }
            
            return true;
        } catch (const std::exception& e) {
            LOG_ERROR("Error saving instance to file: " + std::string(e.what()));
            return false;
        }
    }
};

#endif // DNA_INSTANCE_IO_H 