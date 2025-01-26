#include "../../include/dna/dna_generator.h"
#include "../../include/utils/logging.h"
#include <algorithm>
#include <stdexcept>

DNAGenerator::DNAGenerator(std::unique_ptr<Random> random)
    : m_random(std::move(random)), m_n(0), m_k(0), m_deltaK(0) {
    if (!m_random) {
        m_random = std::make_unique<Random>();
    }
}

void DNAGenerator::setParameters(int n, int k, int deltaK) {
    if (n <= 0) {
        throw std::invalid_argument("DNA length must be positive");
    }
    if (k <= 0 || k > n) {
        throw std::invalid_argument("k-mer length must be positive and not greater than DNA length");
    }
    if (deltaK < 0) {
        throw std::invalid_argument("Delta K must be non-negative");
    }
    m_n = n;
    m_k = k;
    m_deltaK = deltaK;
}

std::string DNAGenerator::generateDNA(int length, bool repAllowed) const {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (length <= 0) {
        throw std::invalid_argument("DNA length must be positive");
    }

    std::string dna;
    dna.reserve(length);

    const char nucleotides[] = {'A', 'C', 'G', 'T'};

    for (int i = 0; i < length; ++i) {
        int index = m_random->getRandomInt(0, 3);
        dna.push_back(nucleotides[index]);
    }

    if (!repAllowed) {
        // Here we could add logic to prevent repetitions if needed
        // For now, we'll just return the DNA as is
        LOG_INFO("Repetition prevention not implemented");
    }

    return dna;
}

bool DNAGenerator::validateParameters() const {
    return m_n > 0 && m_k > 0 && m_k <= m_n && m_deltaK >= 0;
}

DNAInstance DNAGenerator::generateRandomInstance(
    int size,
    int k,
    int lNeg,
    int lPoz,
    int maxErrors,
    bool allowNegative,
    double errorProb) const {
    
    std::lock_guard<std::mutex> lock(m_mutex);

    if (size <= 0 || k <= 0 || k > size) {
        throw std::invalid_argument("Invalid size or k parameters");
    }
    if (lNeg < 0 || lPoz < 0) {
        throw std::invalid_argument("Error counts must be non-negative");
    }
    if (maxErrors < 0) {
        throw std::invalid_argument("Max errors must be non-negative");
    }
    if (errorProb < 0.0 || errorProb > 1.0) {
        throw std::invalid_argument("Error probability must be between 0 and 1");
    }

    DNAInstance instance;
    instance.setN(size);
    instance.setK(k);
    instance.setLNeg(lNeg);
    instance.setLPoz(lPoz);

    // Generate original DNA sequence
    std::string dna = generateDNA(size, true);
    instance.setOriginalDNA(dna);
    instance.setDNA(dna);

    // Generate spectrum
    auto spectrum = generateDNASpectrum(instance);
    instance.setSpectrum(spectrum);

    return instance;
}

bool DNAGenerator::saveToFile(const DNAInstance& instance, const std::string& filename) const {
    try {
        DNAInstanceIO::saveToFile(instance, filename);
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to save instance to file: " + std::string(e.what()));
        return false;
    }
}

DNAInstance DNAGenerator::loadFromFile(const std::string& filename) {
    return DNAInstanceIO::loadFromFile(filename);
}

std::vector<std::string> DNAGenerator::generateDNASpectrum(const DNAInstance& instance) {
    SpectrumGenerator specGen;
    return specGen.generateSpectrum(instance.getDNA(), instance.getK(), instance.getDeltaK());
} 