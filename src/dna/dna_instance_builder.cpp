#include "dna/dna_instance_builder.h"
#include "generator/dna_generator.h"
#include "utils/logging.h"

DNAInstanceBuilder& DNAInstanceBuilder::setN(int value) {
    n = value;
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setK(int value) {
    instance.setK(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setDeltaK(int value) {
    instance.setDeltaK(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLNeg(int value) {
    instance.setLNeg(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLPoz(int value) {
    instance.setLPoz(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setRepAllowed(bool value) {
    instance.setRepAllowed(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setProbablePositive(int value) {
    instance.setProbablePositive(value);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::buildDNA() {
    DNAGenerator generator;
    generator.setParameters(n, instance.getK(), instance.getDeltaK());
    std::string dna = generator.generateDNA(n, instance.isRepAllowed());
    instance.setDNA(dna);
    instance.setN(n);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::buildSpectrum() {
    SpectrumGenerator generator;
    std::vector<std::string> spectrum = generator.generateSpectrum(
        instance.getDNA(), 
        instance.getK(), 
        instance.getDeltaK()
    );
    instance.setSpectrum(spectrum);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::applyError(IErrorIntroductionStrategy* strategy) {
    if (strategy) {
        strategy->introduceErrors(instance);
    }
    return *this;
} 