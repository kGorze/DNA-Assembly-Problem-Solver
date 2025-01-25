#ifndef DNA_INSTANCE_BUILDER_H
#define DNA_INSTANCE_BUILDER_H

#include "dna_instance.h"
#include "generator/dna_generator.h"
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
        n = value;
        return *this;
    }
    
    DNAInstanceBuilder& setK(int value) {
        if (value <= 0) {
            throw std::invalid_argument("K must be positive");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        instance.setK(value);
        return *this;
    }
    
    DNAInstanceBuilder& setDeltaK(int value) {
        if (value < 0) {
            throw std::invalid_argument("DeltaK cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        instance.setDeltaK(value);
        return *this;
    }
    
    DNAInstanceBuilder& setLNeg(int value) {
        if (value < 0) {
            throw std::invalid_argument("LNeg cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        instance.setLNeg(value);
        return *this;
    }
    
    DNAInstanceBuilder& setLPoz(int value) {
        if (value < 0) {
            throw std::invalid_argument("LPoz cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        instance.setLPoz(value);
        return *this;
    }
    
    DNAInstanceBuilder& setRepAllowed(bool value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        instance.setRepAllowed(value);
        return *this;
    }
    
    DNAInstanceBuilder& setProbablePositive(int value) {
        if (value < 0) {
            throw std::invalid_argument("ProbablePositive cannot be negative");
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        instance.setProbablePositive(value);
        return *this;
    }
    
    // Build methods with error handling
    DNAInstanceBuilder& buildDNA() {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            DNAGenerator generator;
            generator.setParameters(n, instance.getK(), instance.getDeltaK());
            std::string dna = generator.generateDNA(n, instance.isRepAllowed());
            if (dna.empty()) {
                throw std::runtime_error("Failed to generate DNA sequence");
            }
            instance.setDNA(dna);
            instance.setN(n);
        } catch (const std::exception& e) {
            LOG_ERROR("Error in buildDNA: " + std::string(e.what()));
            throw;
        }
        return *this;
    }
    
    DNAInstanceBuilder& buildSpectrum() {
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            SpectrumGenerator generator;
            std::vector<std::string> spectrum = generator.generateSpectrum(
                instance.getDNA(), 
                instance.getK(), 
                instance.getDeltaK()
            );
            if (spectrum.empty()) {
                throw std::runtime_error("Failed to generate spectrum");
            }
            instance.setSpectrum(spectrum);
        } catch (const std::exception& e) {
            LOG_ERROR("Error in buildSpectrum: " + std::string(e.what()));
            throw;
        }
        return *this;
    }
    
    DNAInstanceBuilder& applyError(IErrorIntroductionStrategy* strategy) {
        if (!strategy) {
            LOG_WARNING("Null error introduction strategy provided");
            return *this;
        }
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            strategy->introduceErrors(instance);
        } catch (const std::exception& e) {
            LOG_ERROR("Error in applyError: " + std::string(e.what()));
            throw;
        }
        return *this;
    }
    
    // Get the built instance
    DNAInstance getInstance() const { 
        std::lock_guard<std::mutex> lock(m_mutex);
        return instance; 
    }
    
private:
    DNAInstance instance;
    int n = 0;  // DNA length
    ErrorContext errorCtx;
    mutable std::mutex m_mutex;
};

#endif // DNA_INSTANCE_BUILDER_H 