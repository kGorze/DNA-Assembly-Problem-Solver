//
// Created by konrad_guest on 28/12/2024.
// SMART

#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <random>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <mutex>
#include <set>

/* **********************************************
 *          RandomGenerator – singleton
 * **********************************************/
RandomGenerator::RandomGenerator()
{
    std::random_device rd;
    gen.seed(rd() ^ (unsigned long long)std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

RandomGenerator& RandomGenerator::getInstance()
{
    static RandomGenerator instance;
    return instance;
}

std::mt19937& RandomGenerator::get()
{
    return gen;
}

/* **********************************************
 *          DNAGenerator
 * **********************************************/
bool DNAGenerator::validateParameters() const {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_n <= 0) {
        LOG_ERROR("Invalid n value: " + std::to_string(m_n));
        return false;
    }
    if (m_k <= 0) {
        LOG_ERROR("Invalid k value: " + std::to_string(m_k));
        return false;
    }
    if (m_deltaK < 0) {
        LOG_ERROR("Invalid deltaK value: " + std::to_string(m_deltaK));
        return false;
    }
    return true;
}

void DNAGenerator::setParameters(int n, int k, int deltaK) {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (n <= 0) throw std::invalid_argument("n must be positive");
    if (k <= 0) throw std::invalid_argument("k must be positive");
    if (deltaK < 0) throw std::invalid_argument("deltaK cannot be negative");
    
    m_n = n;
    m_k = k;
    m_deltaK = deltaK;
}

std::string DNAGenerator::generateDNA(int length, bool repAllowed) const {
    std::lock_guard<std::mutex> lock(m_mutex);
    try {
        if (length <= 0) {
            throw std::invalid_argument("DNA length must be positive");
        }
        
        if (!m_random) {
            throw std::runtime_error("Random generator not initialized");
        }
        
        const char nucleotides[] = {'A', 'C', 'G', 'T'};
        std::string dna;
        dna.reserve(length);
        
        for (int i = 0; i < length; ++i) {
            int index = (*m_random).getRandomInt(0, 3);
            dna.push_back(nucleotides[index]);
        }
        
        if (!repAllowed) {
            LOG_WARNING("Repetition prevention not implemented");
        }
        
        return dna;
    } catch (const std::exception& e) {
        LOG_ERROR("Error generating DNA: " + std::string(e.what()));
        throw;
    }
}

/* **********************************************
 *          SpectrumGenerator
 * **********************************************/
std::vector<std::string> SpectrumGenerator::generateSpectrum(const std::string &dna, int k, int deltaK)
{
    std::vector<std::string> spectrum;
    if ((int)dna.size() < k) {
        return spectrum;
    }

    auto &rng = RandomGenerator::getInstance().get();
    std::uniform_int_distribution<int> distDelta(0, deltaK);
    std::uniform_int_distribution<int> distSign(0, 1);

    int startPos = 0;
    int currentOligoLength = k;

    // Calculate fixed zone start (last k+2 oligonucleotides)
    int endFixedZoneStart = dna.size() - (k + 2);
    if(endFixedZoneStart < 0) {
        endFixedZoneStart = 0;
    }

    while (startPos + k <= (int)dna.size()) {
        if (startPos <= endFixedZoneStart) {
            if (startPos == 0) {
                // First oligo always has length k
                currentOligoLength = k;
            } else {
                int d = distDelta(rng);
                if (d > 0) {
                    int sign = distSign(rng);
                    // For negative delta, ensure minimum length is k-d
                    if (sign == 0) {
                        currentOligoLength = k - d;
                        // Ensure minimum length is at least k-deltaK
                        currentOligoLength = std::max(currentOligoLength, k - deltaK);
                    } else {
                        currentOligoLength = std::min(k + d, (int)dna.size() - startPos);
                    }
                } else {
                    currentOligoLength = k;
                }
            }
        } else {
            // Last k+2 oligonucleotides always have length k
            currentOligoLength = k;
        }

        if (startPos + currentOligoLength > (int)dna.size()) {
            break;
        }

        spectrum.push_back(dna.substr(startPos, currentOligoLength));
        startPos += 1;
    }

    return spectrum;
}

/* **********************************************
 *      NegativeErrorIntroducer
 * **********************************************/
// void NegativeErrorIntroducer::introduceErrors(DNAInstance& instance) { ... }

/* **********************************************
 *      PositiveErrorIntroducer
 * **********************************************/
// void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) { ... }

DNAInstance DNAGenerator::generateRandomInstance(
    int size,
    int k,
    int lNeg,
    int lPoz,
    int maxErrors,
    bool allowNegative,
    double errorProb
) const {
    if (size <= 0 || k <= 0 || lNeg < 0 || lPoz < 0 || maxErrors < 0 || errorProb < 0.0 || errorProb > 1.0) {
        throw std::invalid_argument("Invalid parameters for random instance generation");
    }

    std::random_device rd;
    int seed = rd();
    return DNAInstance(size, k, lNeg, lPoz, maxErrors, allowNegative, errorProb, seed);
}

bool DNAGenerator::saveToFile(const DNAInstance& instance, const std::string& filename) const {
    if (filename.empty()) {
        throw std::invalid_argument("Filename cannot be empty");
    }
    
    std::lock_guard<std::mutex> lock(m_mutex);
    try {
        return InstanceIO::saveInstance(instance, filename);
    } catch (const std::exception& e) {
        LOG_ERROR("Error saving instance to file: " + std::string(e.what()));
        throw;
    }
}

DNAInstance DNAGenerator::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    DNAInstance instance;
    
    int size;
    std::string sequence;
    
    file >> size;
    file >> sequence;
    
    if (file.fail() || sequence.size() != size_t(size)) {
        LOG_ERROR("Failed to load DNA instance from file");
        return DNAInstance();
    }
    
    instance.setSize(size);
    instance.setTargetSequence(sequence);
    
    LOG_INFO("Successfully loaded DNA instance of size " + std::to_string(instance.getSize()));
    return instance;
}

DNAGenerator::DNAGenerator(std::unique_ptr<Random> random) : m_random(std::move(random)) {}

std::vector<std::string> DNAGenerator::generateDNASpectrum(const DNAInstance& instance) {
    std::vector<std::string> spectrum;
    const std::string& dna = instance.getDNA();
    int k = instance.getK();
    
    if (dna.empty()) {
        throw std::invalid_argument("Invalid DNA or k-mer length: DNA sequence is empty");
    }
    
    for (size_t i = 0; i <= dna.length() - k; ++i) {
        spectrum.push_back(dna.substr(i, k));
    }
    
    return spectrum;
}