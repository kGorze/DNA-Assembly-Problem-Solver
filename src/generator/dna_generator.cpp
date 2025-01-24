//
// Created by konrad_guest on 28/12/2024.
// SMART

#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <random>
#include <fstream>
#include <algorithm>
#include <numeric>

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
std::string DNAGenerator::generateDNA(int n, bool repAllowed)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 3);
    
    const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::string dna;
    dna.reserve(n);
    
    for (int i = 0; i < n; ++i) {
        dna += nucleotides[dis(gen)];
    }
    
    return dna;
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
    auto& spectrum = instance.getSpectrum();
    if (spectrum.empty() || lNeg <= 0) return;
    
    int startIndex = instance.getStartIndex();
    if (startIndex < 0) startIndex = 0;
    
    // Create indices vector
    std::vector<size_t> indices(spectrum.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    // Shuffle indices
    auto& gen = RandomGenerator::getInstance().get();
    std::shuffle(indices.begin(), indices.end(), gen);
    
    // Remove up to lNeg elements
    int toRemove = std::min(lNeg, static_cast<int>(spectrum.size()));
    std::vector<std::string> newSpectrum;
    newSpectrum.reserve(spectrum.size() - toRemove);
    
    for (size_t i = toRemove; i < indices.size(); ++i) {
        newSpectrum.push_back(spectrum[indices[i]]);
    }
    
    // Sort the remaining elements
    std::sort(newSpectrum.begin(), newSpectrum.end());
    instance.setSpectrum(newSpectrum);
}

/* **********************************************
 *      PositiveErrorIntroducer
 * **********************************************/
void PositiveErrorIntroducer::introduceErrors(DNAInstance& instance) {
    if (lPoz <= 0) return;
    
    auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    int deltaK = instance.getDeltaK();
    
    // Choose error introduction method based on probablePositive
    if (instance.getProbablePositive() == 0) {
        // Add random sequences
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> lenDis(-deltaK, deltaK);
        std::uniform_int_distribution<> nucDis(0, 3);
        
        const char nucleotides[] = {'A', 'C', 'G', 'T'};
        
        for (int i = 0; i < lPoz; ++i) {
            int len = k + lenDis(gen);
            std::string newKmer;
            newKmer.reserve(len);
            
            for (int j = 0; j < len; ++j) {
                newKmer += nucleotides[nucDis(gen)];
            }
            
            spectrum.push_back(newKmer);
        }
    } else {
        // Modify existing k-mers
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> posDis(0, k-1);
        std::uniform_int_distribution<> nucDis(0, 3);
        
        const char nucleotides[] = {'A', 'C', 'G', 'T'};
        
        for (int i = 0; i < lPoz && !spectrum.empty(); ++i) {
            std::uniform_int_distribution<> kmerDis(0, spectrum.size()-1);
            int idx = kmerDis(gen);
            std::string modified = spectrum[idx];
            
            // Modify two random positions
            int pos1 = posDis(gen);
            int pos2;
            do {
                pos2 = posDis(gen);
            } while (pos2 == pos1);
            
            modified[pos1] = nucleotides[nucDis(gen)];
            modified[pos2] = nucleotides[nucDis(gen)];
            
            spectrum.push_back(modified);
        }
    }
    
    // Sort the spectrum after adding errors
    std::sort(spectrum.begin(), spectrum.end());
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

bool DNAGenerator::generate() {
    if (!validateParameters()) {
        LOG_ERROR("Invalid parameters for DNA generation");
        return false;
    }
    return true;
}

DNAInstance DNAGenerator::generateRandomInstance(int size) {
    DNAInstance instance;
    
    // Set basic parameters
    instance.setSize(size);
    
    // Generate random DNA sequence
    std::string sequence = generateDNA(size, instance.isRepAllowed());
    instance.setDNA(sequence);
    instance.setTargetSequence(sequence);
    
    return instance;
}

bool DNAGenerator::saveToFile(const DNAInstance& instance, const std::string& filename) {
    std::ofstream file(filename);
    if (!file) return false;
    
    file << instance.getSize() << "\n";
    file << instance.getTargetSequence() << "\n";
    return true;
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