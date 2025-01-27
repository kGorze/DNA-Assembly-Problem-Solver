#include "dna/dna_instance_builder.h"
#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <stdexcept>
#include <random>

DNAInstanceBuilder::DNAInstanceBuilder(std::shared_ptr<DNAGenerator> generator) 
    : m_generator(std::move(generator))
{
    if (!m_generator) {
        LOG_ERROR("DNAInstanceBuilder initialized with null generator");
        throw std::invalid_argument("Generator cannot be null");
    }
    LOG_INFO("DNAInstanceBuilder initialized with generator");
}

bool DNAInstanceBuilder::validateState(const std::string& context) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    bool isValid = true;
    std::string errors;

    // Check basic parameters
    if (m_instance.getN() <= 0) {
        errors += "Invalid N value (" + std::to_string(m_instance.getN()) + "); ";
        isValid = false;
    }
    if (m_instance.getK() <= 0) {
        errors += "Invalid K value (" + std::to_string(m_instance.getK()) + "); ";
        isValid = false;
    }
    if (m_instance.getDeltaK() < 0) {
        errors += "Invalid deltaK value (" + std::to_string(m_instance.getDeltaK()) + "); ";
        isValid = false;
    }

    // Check DNA
    const auto& dna = m_instance.getDNA();
    if (dna.empty()) {
        errors += "DNA is empty; ";
        isValid = false;
    } else if (dna.length() != static_cast<size_t>(m_instance.getN())) {
        errors += "DNA length (" + std::to_string(dna.length()) + 
                 ") doesn't match N (" + std::to_string(m_instance.getN()) + "); ";
        isValid = false;
    }

    // Check spectrum only after buildSpectrum or during final build
    if (context != "buildDNA") {
        const auto& spectrum = m_instance.getSpectrum();
        if (spectrum.empty()) {
            errors += "Spectrum is empty; ";
            isValid = false;
        } else {
            for (const auto& kmer : spectrum) {
                if (kmer.length() != static_cast<size_t>(m_instance.getK())) {
                    errors += "Invalid k-mer length in spectrum: " + std::to_string(kmer.length()) + 
                             " (expected " + std::to_string(m_instance.getK()) + "); ";
                    isValid = false;
                    break;
                }
            }
        }
    }

    if (!isValid) {
        LOG_ERROR("Validation failed in " + context + ": " + errors);
    } else {
        const auto& spectrum = m_instance.getSpectrum();
        LOG_DEBUG("Validation passed in " + context + 
                 " - DNA length: " + std::to_string(dna.length()) + 
                 ", Spectrum size: " + std::to_string(spectrum.size()) + 
                 ", K: " + std::to_string(m_instance.getK()));
    }

    return isValid;
}

DNAInstanceBuilder& DNAInstanceBuilder::setN(int value) {
    m_instance.setN(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setK(int value) {
    m_instance.setK(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setDeltaK(int value) {
    m_instance.setDeltaK(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLNeg(int value) {
    m_instance.setLNeg(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLPoz(int value) {
    m_instance.setLPoz(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setRepAllowed(bool value) {
    m_instance.setRepAllowed(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setProbablePositive(double value) {
    m_instance.setProbablePositive(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setStartIndex(int value) {
    m_instance.setStartIndex(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setDNA(const std::string& dna) {
    m_instance.setDNA(dna);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setSpectrum(const std::vector<std::string>& spectrum) {
    m_instance.setSpectrum(spectrum);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::buildDNA() {
    try {
        LOG_INFO("Building DNA...");
        
        // Validate basic parameters first
        if (m_instance.getN() <= 0) {
            throw std::invalid_argument("DNA length (N) must be positive");
        }
        if (m_instance.getK() <= 0) {
            throw std::invalid_argument("K must be positive");
        }
        if (m_instance.getDeltaK() < 0) {
            throw std::invalid_argument("DeltaK cannot be negative");
        }

        // Set parameters and generate DNA
        m_generator->setParameters(m_instance.getN(), m_instance.getK(), m_instance.getDeltaK());
        std::string dna = m_generator->generateDNA(m_instance.getN(), m_instance.isRepAllowed());
        
        if (dna.empty()) {
            throw std::runtime_error("Generated DNA is empty");
        }
        
        // Set DNA and related fields
        m_instance.setDNA(dna);
        m_instance.setOriginalDNA(dna);
        m_instance.setTargetSequence(dna);
        m_instance.setSize(dna.length());
        
        LOG_INFO("Generated DNA of length " + std::to_string(dna.length()));

        // Generate spectrum immediately after DNA
        LOG_INFO("Generating initial spectrum...");
        auto spectrum = m_generator->generateDNASpectrum(m_instance);
        
        if (!spectrum.empty()) {
            m_instance.setSpectrum(spectrum);
            LOG_INFO("Generated initial spectrum with " + std::to_string(spectrum.size()) + " oligonucleotides");
        }
        
        // Now validate the complete state
        validateState("buildDNA");
        return *this;
    } catch (const std::exception& e) {
        LOG_ERROR("Error building DNA: " + std::string(e.what()));
        throw;
    }
}

DNAInstanceBuilder& DNAInstanceBuilder::buildSpectrum() {
    try {
        LOG_INFO("Building spectrum...");
        validateState("buildSpectrum-pre");
        
        // Generate spectrum using the generator
        LOG_INFO("Generating spectrum...");
        auto spectrum = m_generator->generateDNASpectrum(m_instance);
        
        if (spectrum.empty()) {
            LOG_ERROR("Generated spectrum is empty");
            throw std::runtime_error("Generated spectrum is empty");
        }
        
        // Validate k-mers in spectrum
        for (const auto& kmer : spectrum) {
            if (kmer.length() != static_cast<size_t>(m_instance.getK())) {
                LOG_ERROR("Invalid k-mer length: " + std::to_string(kmer.length()));
                throw std::runtime_error("Invalid k-mer length in spectrum");
            }
        }
        
        m_instance.setSpectrum(spectrum);
        LOG_INFO("Generated spectrum with " + std::to_string(spectrum.size()) + " oligonucleotides");
        LOG_INFO("First k-mer: " + spectrum.front() + ", Last k-mer: " + spectrum.back());
        
        validateState("buildSpectrum-post");
        return *this;
    } catch (const std::exception& e) {
        LOG_ERROR("Error building spectrum: " + std::string(e.what()));
        throw;
    }
}

DNAInstanceBuilder& DNAInstanceBuilder::applyError(IErrorIntroductionStrategy* strategy) {
    if (strategy) {
        try {
            LOG_INFO("Applying error strategy...");
            validateState("applyError-pre");
            
            // Ensure we have a valid spectrum before applying errors
            if (m_instance.getSpectrum().empty()) {
                LOG_INFO("Generating spectrum before applying errors...");
                buildSpectrum();
            }
            
            strategy->introduceErrors(m_instance);
            
            validateState("applyError-post");
            LOG_INFO("Successfully applied error strategy. Final spectrum size: " + 
                    std::to_string(m_instance.getSpectrum().size()));
        } catch (const std::exception& e) {
            LOG_ERROR("Error applying error strategy: " + std::string(e.what()));
            throw;
        }
    }
    return *this;
}

DNAInstance DNAInstanceBuilder::build() const {
    try {
        LOG_INFO("Building final DNA instance...");
        validateState("build-pre");
        
        // Create a new instance
        DNAInstance result;
        
        // Set basic parameters
        result.setN(m_instance.getN());
        result.setK(m_instance.getK());
        result.setDeltaK(m_instance.getDeltaK());
        result.setLNeg(m_instance.getLNeg());
        result.setLPoz(m_instance.getLPoz());
        result.setRepAllowed(m_instance.isRepAllowed());
        result.setProbablePositive(m_instance.getProbablePositive());
        result.setStartIndex(m_instance.getStartIndex());
        result.setSize(m_instance.getSize());
        
        // Set DNA
        result.setDNA(m_instance.getDNA());
        result.setOriginalDNA(m_instance.getOriginalDNA());
        result.setTargetSequence(m_instance.getTargetSequence());
        
        // Set spectrum
        result.setSpectrum(m_instance.getSpectrum());
        
        // Final validation
        const auto& finalSpectrum = result.getSpectrum();
        if (finalSpectrum.empty()) {
            LOG_ERROR("Build failed: final spectrum is empty");
            throw std::runtime_error("Build failed: final spectrum is empty");
        }
        
        LOG_INFO("Successfully built DNA instance");
        LOG_INFO("Final state - DNA length: " + std::to_string(result.getDNA().length()) + 
                ", Spectrum size: " + std::to_string(finalSpectrum.size()) + 
                ", K: " + std::to_string(result.getK()));
        
        return result;
    } catch (const std::exception& e) {
        LOG_ERROR("Error in build: " + std::string(e.what()));
        throw;
    }
} 