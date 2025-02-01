#include "../../include/generator/dna_generator.h"
#include "../../include/utils/logging.h"
#include "../../include/dna/dna_instance_io.h"
#include <algorithm>
#include <stdexcept>
#include <random>

DNAGenerator::DNAGenerator(std::unique_ptr<Random> random) 
    : m_random(std::move(random)), m_n(0), m_k(0), m_deltaK(0) {
    if (!m_random) {
        throw std::invalid_argument("Random generator cannot be null");
    }
}

void DNAGenerator::setParameters(int n, int k, int deltaK) {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    if (n <= 0) throw std::invalid_argument("n must be positive");
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (deltaK < 0) throw std::invalid_argument("deltaK must be non-negative");
    if (k > n) throw std::invalid_argument("k cannot be greater than n");
    
    m_n = n;
    m_k = k;
    m_deltaK = deltaK;
}

std::string DNAGenerator::generateDNA(int length, bool repAllowed) const {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    LOG_INFO("Starting DNA generation with length={}", length);
    
    if (length <= 0) {
        LOG_ERROR("Invalid length: {}", length);
        throw std::invalid_argument("Length must be positive");
    }
    
    static const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::string dna;
    dna.reserve(length);
    
    LOG_INFO("Generating DNA sequence...");
    for (int i = 0; i < length; ++i) {
        if (i % 100 == 0) {
            LOG_INFO("Generated {} nucleotides so far", i);
        }
        int index = m_random->getRandomInt(0, 3);
        dna += nucleotides[index];
    }
    
    if (!repAllowed) {
        LOG_INFO("Repetition prevention not implemented");
    }
    
    LOG_INFO("DNA generation completed, length={}", dna.length());
    return dna;
}

bool DNAGenerator::validateParameters() const {
    return m_n > 0 && m_k > 0 && m_deltaK >= 0 && m_k <= m_n;
}

DNAInstance DNAGenerator::generateRandomInstance(int size, int k, int lNeg, int lPoz) const {
    std::lock_guard<std::recursive_mutex> lock(m_mutex);
    
    LOG_INFO("Starting random instance generation with size={}, k={}, lNeg={}, lPoz={}", size, k, lNeg, lPoz);
    
    if (!validateParameters()) {
        LOG_ERROR("Invalid parameters for DNA generation");
        throw std::invalid_argument("Invalid parameters for DNA generation");
    }

    LOG_INFO("Creating instance with parameters...");
    // Create instance with all required parameters
    DNAInstance instance(size, k, lNeg, lPoz, m_deltaK, false, 0.0, 0);
    
    LOG_INFO("Generating DNA sequence...");
    // Generate and set DNA
    std::string dna = generateDNA(size);
    LOG_INFO("DNA sequence generated, length={}", dna.length());
    
    instance.setDNA(dna);
    instance.setOriginalDNA(dna);
    instance.setTargetSequence(dna);
    
    LOG_INFO("Random instance generation completed");
    return instance;
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