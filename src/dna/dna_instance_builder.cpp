#include "dna/dna_instance_builder.h"
#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <stdexcept>

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

DNAInstanceBuilder& DNAInstanceBuilder::setProbablePositive(int value) {
    m_instance.setProbablePositive(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::buildDNA() {
    try {
        m_generator.setParameters(m_instance.getN(), m_instance.getK(), m_instance.getDeltaK());
        std::string dna = m_generator.generateDNA(m_instance.getN(), m_instance.isRepAllowed());
        m_instance.setDNA(dna);
        return *this;
    } catch (const std::exception& e) {
        LOG_ERROR("Error building DNA: " + std::string(e.what()));
        throw;
    }
}

DNAInstanceBuilder& DNAInstanceBuilder::buildSpectrum() {
    try {
        SpectrumGenerator specGen;
        auto spectrum = specGen.generateSpectrum(
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
    DNAInstance result(std::move(m_instance));
    return result;
} 