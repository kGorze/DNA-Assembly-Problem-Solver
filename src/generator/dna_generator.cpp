//
// Created by konrad_guest on 28/12/2024.
// SMART

#include "generator/dna_generator.h"
#include "utils/logging.h"
#include <random>
#include <fstream>

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
    static const std::string nucleotides = "ACGT";
    auto &rng = RandomGenerator::getInstance().get();
    std::uniform_int_distribution<int> dist(0, 3);

    std::string dna;
    dna.reserve(n);

    // Prosty wariant – generujemy dowolne znaki ACGT.
    // Jeśli repAllowed = false, można np. kontrolować,
    // by nie było 2-3 identycznych znaków pod rząd, itp.
    for(int i = 0; i < n; ++i) {
        dna.push_back(nucleotides[dist(rng)]);
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
 *          InstanceIO
 * **********************************************/
bool InstanceIO::saveInstance(const DNAInstance &instance, const std::string &filename)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        return false;
    }

    out << instance.getN() << " " 
        << instance.getK() << " "
        << instance.getDeltaK() << " "
        << instance.getLNeg() << " "
        << instance.getLPoz() << " "
        << (instance.isRepAllowed() ? 1 : 0) << " "
        << instance.getProbablePositive() << "\n";

    out << instance.getDNA() << "\n";

    out << instance.getSpectrum().size() << "\n";
    for (const auto &frag : instance.getSpectrum()) {
        out << frag << "\n";
    }

    return true;
}

bool InstanceIO::loadInstance(const std::string &filename, DNAInstance &instance)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        return false;
    }

    int n, k, dk, ln, lp, rep, pp;
    in >> n >> k >> dk >> ln >> lp >> rep >> pp;

    instance.setN(n);
    instance.setK(k);
    instance.setDeltaK(dk);
    instance.setLNeg(ln);
    instance.setLPoz(lp);
    instance.setRepAllowed(rep != 0);
    instance.setProbablePositive(pp);

    std::string dna;
    in >> dna;
    instance.setDNA(dna);

    int count;
    in >> count;
    std::vector<std::string> spectrum;
    spectrum.reserve(count);

    for (int i = 0; i < count; i++) {
        std::string frag;
        in >> frag;
        spectrum.push_back(frag);
    }
    instance.setSpectrum(spectrum);

    return true;
}

/* **********************************************
 *          DNAInstanceBuilder
 * **********************************************/
DNAInstanceBuilder& DNAInstanceBuilder::setN(int n)
{
    if(n < 300 || n > 700) {
        n = 400;
    }
    instance.setN(n);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setK(int k)
{
    if(k < 7 || k > 10) {
        k = 8;
    }
    instance.setK(k);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setDeltaK(int dk)
{
    if(dk < 0 || dk > 2) {
        dk = 2;
    }
    instance.setDeltaK(dk);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLNeg(int ln)
{
    if(ln != 0 && ln < 10) {
        ln = 0;
    }
    instance.setLNeg(ln);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLPoz(int lp) {
    if(lp != 0) {
        if(lp < 10) lp = 10;

        if(!instance.getSpectrum().empty()) {
            int maxErrors = instance.getSpectrum().size() / 2;
            if(lp > maxErrors) lp = maxErrors;
        }
    }
    instance.setLPoz(lp);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setRepAllowed(bool rep)
{
    if(instance.getLNeg() > 0) {
        rep = true;
    }
    instance.setRepAllowed(rep);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setProbablePositive(int val)
{
    if(val != 0 && val != 1) {
        val = 0;
    }
    instance.setProbablePositive(val);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::buildDNA()
{
    DNAGenerator dnaGen;
    std::string dna = dnaGen.generateDNA(instance.getN(), instance.isRepAllowed());
    instance.setDNA(dna);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::buildSpectrum()
{
    SpectrumGenerator specGen;
    auto spec = specGen.generateSpectrum(
        instance.getDNA(),
        instance.getK(),
        instance.getDeltaK()
    );
    instance.setSpectrum(spec);

    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::applyError(IErrorIntroductionStrategy* strategy)
{
    errorCtx.setStrategy(strategy);
    errorCtx.execute(instance);
    return *this;
}

DNAInstance DNAInstanceBuilder::getInstance()
{
    return instance;
}

/* **********************************************
 *      NegativeErrorIntroducer
 * **********************************************/
void NegativeErrorIntroducer::introduceErrors(DNAInstance &instance)
{
    auto &spectrum = instance.getSpectrum();
    if(spectrum.empty() || lNeg <= 0) {
        return;
    }
    // Ile fragmentów chcemy usunąć
    int toRemove = std::min((int)spectrum.size() - 1, lNeg); 
    // ^^^ odejmij 1, by nie usuwać startIndex, w razie w.
    
    int startIndex = instance.getStartIndex();
    if (startIndex < 0 || startIndex >= (int)spectrum.size()) {
        // Bezpiecznie: jeżeli startIndex jest out-of-range, wyjdź
        return;
    }

    // Wektor do usunięcia fragmentów (bez startIndex):
    std::vector<int> candidates;
    candidates.reserve(spectrum.size()-1);

    for (int i = 0; i < (int)spectrum.size(); i++) {
        if (i != startIndex) {
            candidates.push_back(i);
        }
    }
    // Tasujemy indeksy kandydatów
    std::shuffle(candidates.begin(), candidates.end(), RandomGenerator::getInstance().get());

    // Usuwamy toRemove pierwszy(ych) z potasowanej listy
    // Pozostałe dajemy do newSpectrum
    std::unordered_set<int> toRemoveSet;
    for (int i = 0; i < toRemove; i++) {
        toRemoveSet.insert(candidates[i]);
    }

    std::vector<std::string> newSpectrum;
    newSpectrum.reserve(spectrum.size() - toRemove);
    // Zostaw startIndex zawsze
    // (możesz go dodać w tej samej kolejności, co był)
    for (int i = 0; i < (int)spectrum.size(); i++) {
        if (toRemoveSet.find(i) == toRemoveSet.end()) {
            newSpectrum.push_back(spectrum[i]);
        }
    }
    instance.setSpectrum(newSpectrum);
}

/* **********************************************
 *      PositiveErrorIntroducer
 * **********************************************/
void PositiveErrorIntroducer::introduceErrors(DNAInstance &instance) {
    auto &spectrum = instance.getSpectrum();
    if(lPoz <= 0 || spectrum.empty()) {
        return;
    }

    int k = instance.getK();
    int deltaK = instance.getDeltaK();
    int probablePositive = instance.getProbablePositive();
    
    auto &rng = RandomGenerator::getInstance().get();
    std::uniform_int_distribution<int> distDelta(0, deltaK);
    std::uniform_int_distribution<int> distSign(0, 1);
    std::uniform_int_distribution<int> nucDist(0, 3);
    
    static const std::string nucleotides = "ACGT";
    std::unordered_set<std::string> uniqueSpectrum(spectrum.begin(), spectrum.end());
    
    auto generateRandomKmer = [&](int length) {
        std::string kmer;
        kmer.reserve(length);
        for(int i=0; i<length; i++) {
            kmer.push_back(nucleotides[nucDist(rng)]);
        }
        return kmer;
    };
    
    int errorsAdded = 0;
    int maxAttempts = lPoz * 10;
    int attempts = 0;
    
    while(errorsAdded < lPoz && attempts < maxAttempts) {
        attempts++;
        
        if(probablePositive == 0) {
            int d = distDelta(rng);
            int length = k;
            if(d > 0) {
                length = distSign(rng) ? k + d : k - d;
                length = std::max(1, std::min(length, k + deltaK));
            }
            
            std::string newKmer = generateRandomKmer(length);
            if(uniqueSpectrum.insert(newKmer).second) {
                spectrum.push_back(newKmer);
                errorsAdded++;
            }
        } else {
            std::uniform_int_distribution<int> indexDist(0, spectrum.size() - 1);
            int idx = indexDist(rng);
            std::string original = spectrum[idx];
            
            std::string modified = original;
            int posToChange = (distSign(rng) == 0) ? 
                (int)modified.length() - 1 : 
                (int)modified.length() / 2;
                
            char originalNuc = modified[posToChange];
            char newNuc;
            do {
                newNuc = nucleotides[nucDist(rng)];
            } while(newNuc == originalNuc);
            
            modified[posToChange] = newNuc;
            
            if(uniqueSpectrum.insert(modified).second) {
                spectrum.push_back(modified);
                errorsAdded++;
            }
        }
    }
    
    if(errorsAdded < lPoz) {
        LOG_WARNING("Warning: Could only add " + std::to_string(errorsAdded) + 
                     " positive errors out of " + std::to_string(lPoz) + " requested");
    }
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
    LOG_INFO("Generating new DNA instance");
    LOG_DEBUG("Parameters: N=" + std::to_string(n) + 
              ", K=" + std::to_string(k) + 
              ", deltaK=" + std::to_string(deltaK));
    
    if (!validateParameters()) {
        return false;
    }
    // Add implementation here
    return true;  // Add return statement
}

DNAInstance DNAGenerator::generateRandomInstance(int size) {
    LOG_INFO("Generating random DNA instance of size " + std::to_string(size));
    
    if (size <= 0) {
        LOG_ERROR("Invalid size for DNA instance: " + std::to_string(size));
        return DNAInstance();
    }
    
    static std::mt19937 rng(std::random_device{}());
    static const std::string nucleotides = "ACGT";
    std::uniform_int_distribution<int> dist(0, nucleotides.size() - 1);
    
    DNAInstance instance;
    instance.setSize(size);
    std::string sequence;
    sequence.reserve(size);
    
    for (int i = 0; i < size; i++) {
        sequence.push_back(nucleotides[dist(rng)]);
    }
    instance.setTargetSequence(sequence);
    
    LOG_INFO("Successfully generated DNA instance");
    DEBUG_LOG("First 10 nucleotides: " + 
              instance.getTargetSequence().substr(0, std::min(10, size)));
    
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