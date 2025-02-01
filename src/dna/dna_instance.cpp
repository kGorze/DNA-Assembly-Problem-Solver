#include "../../include/dna/dna_instance.h"
#include "../../include/utils/logging.h"
#include <random>
#include <algorithm>
#include <stdexcept>
#include <sstream>

DNAInstance::DNAInstance(int n, int k, int lNeg, int lPoz, int maxErrors, bool allowNegative, double errorProb, int seed)
    : n(n), k(k), deltaK(maxErrors), lNeg(lNeg), lPoz(lPoz), repAllowed(allowNegative), probablePositive(errorProb), size(n) {
    
    LOG_INFO("Creating DNAInstance with parameters: n={}, k={}, lNeg={}, lPoz={}, maxErrors={}, allowNegative={}, errorProb={}, seed={}", 
             n, k, lNeg, lPoz, maxErrors, allowNegative, errorProb, seed);
    
    if (n <= 0 || k <= 0 || k > n) {
        LOG_ERROR("Invalid parameters: n={}, k={}", n, k);
        throw std::invalid_argument("Invalid DNA or k-mer length");
    }

    LOG_INFO("Creating random generator with seed {}", seed);
    m_random = std::make_unique<Random>(seed);
    
    LOG_INFO("Generating random DNA sequence of length {}", n);
    m_dna = generateRandomDNA(n, *m_random);
    LOG_INFO("Generated DNA: {}", m_dna);
    
    m_originalDNA = m_dna;
    targetSequence = m_dna;  // Initialize target sequence with DNA
    
    LOG_INFO("Generating spectrum...");
    generateSpectrum();
    LOG_INFO("DNAInstance creation completed");
}

void DNAInstance::generateSpectrum() {
    LOG_INFO("Starting spectrum generation...");
    if (k <= 0) {
        LOG_ERROR("Invalid k-mer length: {}", k);
        throw std::invalid_argument("k-mer length must be positive");
    }

    if (static_cast<size_t>(k) > m_dna.length()) {
        LOG_ERROR("k-mer length {} exceeds DNA length {}", k, m_dna.length());
        throw std::invalid_argument("k-mer length exceeds DNA length");
    }

    LOG_INFO("Clearing existing spectrum...");
    m_spectrum.clear();
    m_spectrum.reserve(m_dna.length() - k + 1);

    LOG_INFO("Generating k-mers...");
    // Generate k-mers
    for (size_t i = 0; i + k <= m_dna.length(); ++i) {
        m_spectrum.push_back(m_dna.substr(i, k));
        if (i % 1000 == 0) {
            LOG_INFO("Generated {} k-mers so far...", i);
        }
    }

    LOG_INFO("Sorting spectrum...");
    // Sort spectrum for consistent ordering
    std::sort(m_spectrum.begin(), m_spectrum.end());

    // Remove duplicates if repetitions are not allowed
    if (!repAllowed) {
        LOG_INFO("Removing duplicates...");
        auto last = std::unique(m_spectrum.begin(), m_spectrum.end());
        m_spectrum.erase(last, m_spectrum.end());
    }

    if (!validateSpectrum()) {
        LOG_ERROR("Generated spectrum is invalid");
        throw std::runtime_error("Generated spectrum is invalid");
    }

    LOG_INFO("Spectrum generation completed. Generated {} k-mers", m_spectrum.size());
}

std::string DNAInstance::generateRandomDNA(int length, Random& random) const {
    LOG_INFO("Starting generateRandomDNA with length={}", length);
    
    static const char nucleotides[] = "ACGT";
    static const int numNucleotides = 4;

    std::string dna;
    dna.reserve(length);

    LOG_INFO("Generating DNA sequence...");
    for (int i = 0; i < length; ++i) {
        int index = random.getRandomInt(0, numNucleotides - 1);
        LOG_INFO("Generated index {} for position {}", index, i);
        dna += nucleotides[index];
    }
    
    LOG_INFO("Generated DNA sequence: {}", dna);
    return dna;
}

bool DNAInstance::validateSpectrum() const {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_spectrum.empty()) {
        LOG_ERROR("Spectrum validation failed: spectrum is empty");
        return false;
    }
    
    // Check that all k-mers have the same length
    const size_t kmerLength = m_spectrum[0].length();
    for (const auto& kmer : m_spectrum) {
        if (kmer.length() != kmerLength) {
            LOG_ERROR("Spectrum validation failed: inconsistent k-mer lengths");
            return false;
        }
        // Check that k-mer only contains valid DNA bases
        if (kmer.find_first_not_of("ACGT") != std::string::npos) {
            LOG_ERROR("Spectrum validation failed: invalid DNA bases in k-mer: " + kmer);
            return false;
        }
    }
    
    LOG_DEBUG("Spectrum validation passed: " + std::to_string(m_spectrum.size()) + " k-mers of length " + std::to_string(kmerLength));
    return true;
}

int DNAInstance::findStartVertexIndex(const DNAInstance& instance) {
    const auto& spectrum = instance.getSpectrum();
    if (spectrum.empty()) {
        LOG_ERROR("Cannot find start vertex: empty spectrum");
        return -1;
    }

    // Find the lexicographically smallest k-mer
    auto minKmer = std::min_element(spectrum.begin(), spectrum.end());
    if (minKmer == spectrum.end()) {
        LOG_ERROR("Failed to find minimum k-mer");
        return -1;
    }

    return static_cast<int>(std::distance(spectrum.begin(), minKmer));
}

void DNAInstance::setDNA(const std::string& dna) {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (dna.empty()) {
        LOG_ERROR("Cannot set empty DNA sequence");
        throw std::invalid_argument("DNA sequence cannot be empty");
    }
    
    m_dna = dna;
    size = dna.length();
    
    if (targetSequence.empty()) {
        targetSequence = dna;  // Set target sequence if not already set
    }
    
    if (m_originalDNA.empty()) {
        m_originalDNA = dna;  // Set original DNA if not already set
    }
}