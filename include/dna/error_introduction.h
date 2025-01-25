#pragma once

#include "utils/logging.h"
#include "dna/dna_instance.h"
#include <memory>
#include <random>
#include <mutex>

// Interface for error introduction strategies
class IErrorIntroductionStrategy {
public:
    virtual ~IErrorIntroductionStrategy() = default;
    virtual void introduceErrors(DNAInstance& instance) = 0;
};

// Context for error introduction
struct ErrorContext {
    int lNeg = 0;
    int lPoz = 0;
    bool repAllowed = false;
    int probablePositive = 0;
};

// Base class for error introducers
class BaseErrorIntroducer : public IErrorIntroductionStrategy {
protected:
    mutable std::mutex m_mutex;
    
    // Thread-safe random number generation
    int getRandomNumber(int min, int max) const {
        std::lock_guard<std::mutex> lock(m_mutex);
        static thread_local std::random_device rd;
        static thread_local std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(min, max);
        return dis(gen);
    }
    
    // Validate instance before introducing errors
    bool validateInstance(const DNAInstance& instance) const {
        if (instance.getDNA().empty()) {
            LOG_ERROR("Cannot introduce errors: DNA sequence is empty");
            return false;
        }
        if (instance.getSpectrum().empty()) {
            LOG_ERROR("Cannot introduce errors: spectrum is empty");
            return false;
        }
        return true;
    }
};

// Strategy for introducing negative errors (removing k-mers)
class NegativeErrorIntroducer : public BaseErrorIntroducer {
private:
    int m_lNeg;
    mutable std::mt19937 m_rng;

public:
    explicit NegativeErrorIntroducer(int lNeg) : m_lNeg(lNeg), m_rng(std::random_device{}()) {}
    
    void introduceErrors(DNAInstance& instance) override {
        if (m_lNeg <= 0) return;
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            if (!validateInstance(instance)) {
                throw std::runtime_error("Invalid instance for negative error introduction");
            }
            
            auto spectrum = instance.getSpectrum();
            if (spectrum.size() <= static_cast<size_t>(m_lNeg)) {
                LOG_WARNING("Cannot remove " + std::to_string(m_lNeg) + 
                          " k-mers from spectrum of size " + std::to_string(spectrum.size()));
                return;
            }
            
            // Remove random k-mers
            for (int i = 0; i < m_lNeg; ++i) {
                int idx = getRandomNumber(0, spectrum.size() - 1);
                spectrum.erase(spectrum.begin() + idx);
            }
            
            instance.setSpectrum(spectrum);
        } catch (const std::exception& e) {
            LOG_ERROR("Error in negative error introduction: " + std::string(e.what()));
            throw;
        }
    }
};

// Strategy for introducing positive errors (adding incorrect k-mers)
class PositiveErrorIntroducer : public BaseErrorIntroducer {
private:
    int m_lPoz;
    int m_k;
    mutable std::mt19937 m_rng;

public:
    PositiveErrorIntroducer(int lPoz, int k) : m_lPoz(lPoz), m_k(k), m_rng(std::random_device{}()) {}
    
    void introduceErrors(DNAInstance& instance) override {
        if (m_lPoz <= 0) return;
        
        std::lock_guard<std::mutex> lock(m_mutex);
        try {
            if (!validateInstance(instance)) {
                throw std::runtime_error("Invalid instance for positive error introduction");
            }
            
            auto spectrum = instance.getSpectrum();
            const int deltaK = instance.getDeltaK();
            
            // Add random k-mers
            for (int i = 0; i < m_lPoz; ++i) {
                int len = m_k + getRandomNumber(-deltaK, deltaK);
                std::string newKmer = generateRandomKmer(len);
                spectrum.push_back(newKmer);
            }
            
            // Sort the spectrum after adding errors
            std::sort(spectrum.begin(), spectrum.end());
            instance.setSpectrum(spectrum);
        } catch (const std::exception& e) {
            LOG_ERROR("Error in positive error introduction: " + std::string(e.what()));
            throw;
        }
    }
    
    std::string generateRandomKmer(int length) const {
        static const char nucleotides[] = {'A', 'C', 'G', 'T'};
        std::string kmer;
        kmer.reserve(length);
        
        for (int i = 0; i < length; ++i) {
            kmer.push_back(nucleotides[getRandomNumber(0, 3)]);
        }
        
        return kmer;
    }
};

// Factory for creating error introducers
class ErrorIntroducerFactory {
public:
    static std::unique_ptr<IErrorIntroductionStrategy> createNegativeErrorIntroducer(int lNeg) {
        return std::make_unique<NegativeErrorIntroducer>(lNeg);
    }
    
    static std::unique_ptr<IErrorIntroductionStrategy> createPositiveErrorIntroducer(int lPoz, int k) {
        return std::make_unique<PositiveErrorIntroducer>(lPoz, k);
    }
}; 