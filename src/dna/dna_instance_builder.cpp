#include "dna/dna_instance_builder.h"
#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <stdexcept>
#include <random>

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
        // Validate parameters
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
        
        m_instance.setDNA(dna);
        m_instance.setOriginalDNA(dna);
        m_instance.setTargetSequence(dna);
        m_instance.setSize(dna.length());
        
        LOG_INFO("Successfully generated DNA of length " + std::to_string(dna.length()));
        return *this;
    } catch (const std::exception& e) {
        LOG_ERROR("Error building DNA: " + std::string(e.what()));
        throw;
    }
}

DNAInstanceBuilder& DNAInstanceBuilder::buildSpectrum() {
    try {
        // Validate DNA exists and has valid length
        const std::string& dna = m_instance.getDNA();
        if (dna.empty()) {
            LOG_ERROR("Cannot generate spectrum: DNA is empty");
            throw std::runtime_error("Cannot generate spectrum: DNA is empty");
        }
        LOG_INFO("DNA length for spectrum generation: " + std::to_string(dna.length()));
        
        // Validate parameters
        int k = m_instance.getK();
        int deltaK = m_instance.getDeltaK();
        if (k <= 0) {
            LOG_ERROR("Invalid k value: " + std::to_string(k));
            throw std::invalid_argument("K must be positive");
        }
        if (deltaK < 0) {
            LOG_ERROR("Invalid deltaK value: " + std::to_string(deltaK));
            throw std::invalid_argument("DeltaK cannot be negative");
        }
        if (k > static_cast<int>(dna.length())) {
            LOG_ERROR("k value (" + std::to_string(k) + ") larger than DNA length (" + std::to_string(dna.length()) + ")");
            throw std::invalid_argument("K cannot be larger than DNA length");
        }
        LOG_INFO("Using k=" + std::to_string(k) + ", deltaK=" + std::to_string(deltaK));
        
        // Create spectrum generator and generate spectrum
        SpectrumGenerator specGen;
        LOG_INFO("Generating spectrum...");
        auto spectrum = specGen.generateSpectrum(dna, k, deltaK);
        
        if (spectrum.empty()) {
            LOG_ERROR("Generated spectrum is empty");
            throw std::runtime_error("Generated spectrum is empty");
        }
        
        // Validate k-mers in spectrum
        for (const auto& kmer : spectrum) {
            if (kmer.length() < k - deltaK || kmer.length() > k + deltaK) {
                LOG_ERROR("Invalid k-mer length in spectrum: " + std::to_string(kmer.length()));
                throw std::runtime_error("Invalid k-mer length in spectrum");
            }
        }
        
        m_instance.setSpectrum(spectrum);
        LOG_INFO("Successfully generated spectrum with " + std::to_string(spectrum.size()) + " oligonucleotides");
        LOG_INFO("First k-mer: " + spectrum.front() + ", Last k-mer: " + spectrum.back());
        return *this;
    } catch (const std::exception& e) {
        LOG_ERROR("Error building spectrum: " + std::string(e.what()));
        throw;
    }
}

DNAInstanceBuilder& DNAInstanceBuilder::applyError(IErrorIntroductionStrategy* strategy) {
    if (strategy) {
        try {
            strategy->introduceErrors(m_instance);
        } catch (const std::exception& e) {
            LOG_ERROR("Error applying error strategy: " + std::string(e.what()));
            throw;
        }
    }
    return *this;
}

DNAInstance DNAInstanceBuilder::build() const {
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
    result.setTargetSequence(m_instance.getTargetSequence());
    
    // Set DNA and spectrum
    result.setOriginalDNA(m_instance.getOriginalDNA());
    result.setDNA(m_instance.getDNA());
    result.setSpectrum(m_instance.getSpectrum());
    
    return result;
} 