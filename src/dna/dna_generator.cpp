#include "../../include/generator/dna_generator.h"
#include "../../include/utils/logging.h"
#include "../../include/dna/dna_instance_io.h"
#include <algorithm>
#include <stdexcept>
#include <random>

DNAGenerator::DNAGenerator(std::unique_ptr<Random> random) 
    : m_random(random ? std::move(random) : std::make_unique<Random>()) {}

void DNAGenerator::setParameters(int n, int k, int deltaK) {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (n <= 0) throw std::invalid_argument("n must be positive");
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (deltaK < 0) throw std::invalid_argument("deltaK must be non-negative");
    if (k > n) throw std::invalid_argument("k cannot be greater than n");
    
    m_n = n;
    m_k = k;
    m_deltaK = deltaK;
}

std::string DNAGenerator::generateDNA(int length, bool repAllowed) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (length <= 0) throw std::invalid_argument("Length must be positive");
    
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::string dna;
    dna.reserve(length);
    
    for (int i = 0; i < length; ++i) {
        int index = m_random->getRandomInt(0, 3);
        dna += nucleotides[index];
    }
    
    if (!repAllowed) {
        // Here you could add logic to prevent repetitions if needed
        // For now, we'll just return the sequence as is
        LOG_INFO("Repetition prevention not implemented");
    }
    
    return dna;
}

bool DNAGenerator::validateParameters() const {
    return m_n > 0 && m_k > 0 && m_deltaK >= 0 && m_k <= m_n;
}

DNAInstance DNAGenerator::generateRandomInstance(int size, int k, int lNeg, int lPoz) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    
    if (!validateParameters()) {
        throw std::invalid_argument("Invalid parameters for DNA generation");
    }

    auto dna = generateDNA(size);
    // Create instance with all required parameters
    return DNAInstance(size, k, lNeg, lPoz, 0, false, 0.0, 0);
}

bool DNAGenerator::saveToFile(const DNAInstance& instance, const std::string& filename) const {
    try {
        return InstanceIO::saveInstance(instance, filename);
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to save instance to file: " + std::string(e.what()));
        return false;
    }
}

DNAInstance DNAGenerator::loadFromFile(const std::string& filename) {
    return InstanceIO::loadInstance(filename);
}

std::vector<std::string> DNAGenerator::generateDNASpectrum(const DNAInstance& instance) {
    SpectrumGenerator specGen;
    return specGen.generateSpectrum(instance.getDNA(), instance.getK(), instance.getDeltaK());
} 