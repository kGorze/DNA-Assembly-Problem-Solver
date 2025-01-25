#ifndef DNA_INSTANCE_BUILDER_H
#define DNA_INSTANCE_BUILDER_H

#include "dna_instance.h"
#include "generator/dna_generator.h"
#include "error_introduction.h"
#include <string>
#include <vector>
#include <mutex>
#include <memory>
#include <stdexcept>

class DNAInstanceBuilder {
public:
    DNAInstanceBuilder() = default;
    
    // Builder methods with validation
    DNAInstanceBuilder& setN(int value) {
        if (value <= 0) {
            throw std::invalid_argument("N must be positive");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        m_n = value;
        return *this;
    }
    
    DNAInstanceBuilder& setK(int value) {
        if (value <= 0) {
            throw std::invalid_argument("K must be positive");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        m_instance.setK(value);
        return *this;
    }
    
    DNAInstanceBuilder& setDeltaK(int value) {
        if (value < 0) {
            throw std::invalid_argument("DeltaK cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        m_instance.setDeltaK(value);
        return *this;
    }
    
    DNAInstanceBuilder& setLNeg(int value) {
        if (value < 0) {
            throw std::invalid_argument("LNeg cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        m_instance.setLNeg(value);
        return *this;
    }
    
    DNAInstanceBuilder& setLPoz(int value) {
        if (value < 0) {
            throw std::invalid_argument("LPoz cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        m_instance.setLPoz(value);
        return *this;
    }
    
    DNAInstanceBuilder& setRepAllowed(bool value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_instance.setRepAllowed(value);
        return *this;
    }
    
    DNAInstanceBuilder& setProbablePositive(int value) {
        if (value < 0) {
            throw std::invalid_argument("ProbablePositive cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        m_instance.setProbablePositive(value);
        return *this;
    }
    
    DNAInstanceBuilder& buildDNA() {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            DNAGenerator generator;
            generator.setParameters(m_n, m_instance.getK(), m_instance.getDeltaK());
            std::string dna = generator.generateDNA(m_n, m_instance.isRepAllowed());
            m_instance.setDNA(dna);
            m_instance.setN(m_n);
            return *this;
        } catch (const std::exception& e) {
            LOG_ERROR("Error building DNA: " + std::string(e.what()));
            throw;
        }
    }
    
    DNAInstanceBuilder& buildSpectrum() {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            SpectrumGenerator generator;
            std::vector<std::string> spectrum = generator.generateSpectrum(
                m_instance.getDNA(), 
                m_instance.getK(), 
                m_instance.getDeltaK()
            );
            m_instance.setSpectrum(spectrum);
            return *this;
        } catch (const std::exception& e) {
            LOG_ERROR("Error building spectrum: " + std::string(e.what()));
            throw;
        }
    }
    
    DNAInstanceBuilder& introduceErrors() {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            // Find start index before introducing errors
            std::string startFrag = m_instance.getDNA().substr(0, m_instance.getK());
            const auto& spectrum = m_instance.getSpectrum();
            
            int startIdx = -1;
            for (int i = 0; i < static_cast<int>(spectrum.size()); i++) {
                if (spectrum[i] == startFrag) {
                    startIdx = i;
                    break;
                }
            }
            m_instance.setStartIndex(startIdx);
            
            // Introduce negative errors first
            if (m_instance.getLNeg() > 0) {
                auto negErr = ErrorIntroducerFactory::createNegativeErrorIntroducer(m_instance.getLNeg());
                negErr->introduceErrors(m_instance);
            }
            
            // Then introduce positive errors
            if (m_instance.getLPoz() > 0) {
                auto posErr = ErrorIntroducerFactory::createPositiveErrorIntroducer(m_instance.getLPoz());
                posErr->introduceErrors(m_instance);
            }
            
            return *this;
        } catch (const std::exception& e) {
            LOG_ERROR("Error introducing errors: " + std::string(e.what()));
            throw;
        }
    }
    
    const DNAInstance& getInstance() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_instance;
    }
    
private:
    mutable std::mutex m_mutex;
    DNAInstance m_instance;
    int m_n = 0;
};

#endif // DNA_INSTANCE_BUILDER_H 