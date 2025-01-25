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
        
        static thread_local std::random_device rd;
        static thread_local std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 3);
        
        const char nucleotides[] = {'A', 'C', 'G', 'T'};
        std::string dna;
        dna.reserve(length);
        
        for (int i = 0; i < length; ++i) {
            dna.push_back(nucleotides[dis(gen)]);
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
void NegativeErrorIntroducer::introduceErrors(DNAInstance& instance) {
    auto spectrum = instance.getSpectrum();
    if (spectrum.empty() || m_lNeg <= 0) return;

    std::uniform_int_distribution<int> dist(0, spectrum.size() - 1);
    std::vector<int> indices;
    indices.reserve(spectrum.size());
    for (size_t i = 0; i < spectrum.size(); ++i) {
        indices.push_back(i);
    }

    std::shuffle(indices.begin(), indices.end(), m_rng);

    int toRemove = std::min(m_lNeg, static_cast<int>(spectrum.size()));
    indices.resize(toRemove);

    std::sort(indices.begin(), indices.end(), std::greater<int>());
    for (int idx : indices) {
        spectrum.erase(spectrum.begin() + idx);
    }

    instance.setSpectrum(spectrum);
}

/* **********************************************
 *      PositiveErrorIntroducer
 * **********************************************/
void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (m_lPoz <= 0) return;

    auto spectrum = instance.getSpectrum();
    auto k = instance.getK();
    auto deltaK = instance.getDeltaK();

    std::uniform_int_distribution<int> lengthDist(k - deltaK, k + deltaK);
    std::uniform_int_distribution<int> baseDist(0, 3);

    if (!instance.isRepAllowed()) {
        std::set<std::string> uniqueOligos;
        for (const auto& oligo : spectrum) {
            uniqueOligos.insert(oligo);
        }

        for (int i = 0; i < m_lPoz && !spectrum.empty(); ++i) {
            int length = lengthDist(m_rng);
            std::string oligo;
            oligo.reserve(length);
            bool found = false;
            int attempts = 0;
            const int maxAttempts = 100;

            while (!found && attempts < maxAttempts) {
                oligo.clear();
                for (int j = 0; j < length; ++j) {
                    oligo += "ACGT"[baseDist(m_rng)];
                }
                if (uniqueOligos.find(oligo) == uniqueOligos.end()) {
                    found = true;
                    uniqueOligos.insert(oligo);
                    spectrum.push_back(oligo);
                }
                ++attempts;
            }
        }
    } else {
        for (int i = 0; i < m_lPoz; ++i) {
            int length = lengthDist(m_rng);
            std::string oligo;
            oligo.reserve(length);
            for (int j = 0; j < length; ++j) {
                oligo += "ACGT"[baseDist(m_rng)];
            }
            spectrum.push_back(oligo);
        }
    }

    instance.setSpectrum(spectrum);
}

int DNAInstance::findStartVertexIndex(const DNAInstance& instance) {
    const std::string& startFragment = instance.getDNA().substr(0, instance.getK());
    const auto& spectrum = instance.getSpectrum();

    // Znajdź indeks fragmentu startowego w spektrum
    auto it = std::find(spectrum.begin(), spectrum.end(), startFragment);
    
    if (it != spectrum.end()) {
        return std::distance(spectrum.begin(), it);
    }

    // Jeśli nie znaleziono, zwróć -1 lub obsłuż błąd
    return -1;
}

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