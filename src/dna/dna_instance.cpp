#include "dna/dna_instance.h"
#include "utils/logging.h"
#include <algorithm>
#include <random>
#include <chrono>

std::string DNAInstance::generateRandomDNA(int length, Random& random) const {
    std::string dna;
    dna.reserve(length);
    const std::string nucleotides = "ACGT";
    
    for (int i = 0; i < length; ++i) {
        dna += nucleotides[random.getGenerator()() % 4];
    }
    
    return dna;
}

void DNAInstance::generateSpectrum() {
    if (m_originalDNA.empty() || k <= 0) {
        LOG_WARNING("Cannot generate spectrum: DNA is empty or k <= 0");
        return;
    }

    m_spectrum.clear();
    
    for (size_t i = 0; i + k <= m_originalDNA.length(); ++i) {
        m_spectrum.push_back(m_originalDNA.substr(i, k));
    }
    
    if (!repAllowed) {
        std::sort(m_spectrum.begin(), m_spectrum.end());
        m_spectrum.erase(std::unique(m_spectrum.begin(), m_spectrum.end()), m_spectrum.end());
    }
}

int DNAInstance::findStartVertexIndex(const DNAInstance& instance) {
    if (instance.getSpectrum().empty() || instance.getDNA().empty() || instance.getK() <= 0) {
        return -1;
    }

    std::string startFrag = instance.getDNA().substr(0, instance.getK());
    const auto& spectrum = instance.getSpectrum();
    
    for (int i = 0; i < static_cast<int>(spectrum.size()); i++) {
        if (spectrum[i] == startFrag) {
            return i;
        }
    }
    
    return -1;
}

DNAInstance::DNAInstance(int n_val, int k_val, int lNeg_val, int lPoz_val, int maxErrors, bool allowNegative, double errorProb, int seed) {
    n = n_val;
    k = k_val;
    lNeg = lNeg_val;
    lPoz = lPoz_val;
    repAllowed = allowNegative;
    probablePositive = errorProb;
    m_random = std::make_unique<Random>(seed);
    
    // Generate random DNA sequence
    const std::string nucleotides = "ACGT";
    std::uniform_int_distribution<> dist(0, 3);
    
    m_originalDNA.reserve(n);
    for (int i = 0; i < n; ++i) {
        m_originalDNA += nucleotides[dist(m_random->getGenerator())];
    }
    
    // Generate spectrum
    generateSpectrum();
} 