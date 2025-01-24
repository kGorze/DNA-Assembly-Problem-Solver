#ifndef DNA_INSTANCE_H
#define DNA_INSTANCE_H

#include <string>
#include <vector>

class DNAInstance {
public:
    DNAInstance() = default;
    
    // Getters for instance parameters
    int getN() const { return n; }
    int getK() const { return k; }
    int getDeltaK() const { return deltaK; }
    int getLNeg() const { return lNeg; }
    int getLPoz() const { return lPoz; }
    bool isRepAllowed() const { return repAllowed; }
    int getProbablePositive() const { return probablePositive; }
    int getStartIndex() const { return startIndex; }
    int getSize() const { return size; }
    const std::string& getTargetSequence() const { return targetSequence; }
    
    // Getters for DNA and spectrum
    const std::string& getDNA() const { return dna; }
    const std::vector<std::string>& getSpectrum() const { return spectrum; }
    std::vector<std::string>& getSpectrum() { return spectrum; }
    
    // Setters
    void setN(int value) { n = value; }
    void setK(int value) { k = value; }
    void setDeltaK(int value) { deltaK = value; }
    void setLNeg(int value) { lNeg = value; }
    void setLPoz(int value) { lPoz = value; }
    void setRepAllowed(bool value) { repAllowed = value; }
    void setProbablePositive(int value) { probablePositive = value; }
    void setDNA(const std::string& value) { dna = value; }
    void setSpectrum(const std::vector<std::string>& value) { spectrum = value; }
    void setStartIndex(int value) { startIndex = value; }
    void setSize(int value) { size = value; }
    void setTargetSequence(const std::string& value) { targetSequence = value; }

    // Additional functionality
    int findStartVertexIndex(const DNAInstance& instance);

private:
    int n = 0;  // DNA length
    int k = 0;
    int deltaK = 0;
    int lNeg = 0;
    int lPoz = 0;
    bool repAllowed = false;
    int probablePositive = 0;
    int startIndex = -1;
    int size = 0;
    
    std::string dna;
    std::string targetSequence;
    std::vector<std::string> spectrum;
};

#endif // DNA_INSTANCE_H 