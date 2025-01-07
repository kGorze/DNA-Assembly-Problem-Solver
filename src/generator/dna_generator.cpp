//
// Created by konrad_guest on 28/12/2024.
//

#include "generator/dna_generator.h"


bool InstanceIO::saveInstance(const DNAInstance &instance, const std::string &filename)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        return false;
    }

    // Zapisujemy: n, k
    out << instance.getN() << " " << instance.getK() << "\n";

    // Zapisujemy DNA
    out << instance.getDNA() << "\n";

    // Zapisujemy rozmiar spektrum i poszczególne k-mery
    out << instance.getSpectrum().size() << "\n";
    for (const auto &frag : instance.getSpectrum()) {
        out << frag << "\n";
    }

    // Zamykanie pliku przez destruktor ofstream
    return true;
}

bool InstanceIO::loadInstance(const std::string &filename, DNAInstance &instance)
{
    std::ifstream in(filename);
    if (!in.is_open()) {
        return false;
    }

    int n, k;
    in >> n >> k;
    instance.setN(n);
    instance.setK(k);

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

std::string DNAGenerator::generateDNA(int n) {
    // Zakładamy, że mamy 4 nukleotydy: A, C, G, T
    static const std::string nucleotides = "ACGT";

    // Pobranie globalnego generatora pseudolosowego
    auto &rng = RandomGenerator::getInstance().get();
    // Rozkład 0..3
    std::uniform_int_distribution<int> dist(0, 3);

    std::string dna;
    dna.reserve(n);

    for(int i = 0; i < n; ++i) {
        // Losowy indeks i wybieramy literę z "ACGT"
        int idx = dist(rng);
        dna.push_back(nucleotides[idx]);
    }

    return dna;
}

std::vector<std::string> SpectrumGenerator::generateSpectrum(const std::string &dna, int k) {
    std::vector<std::string> spectrum;
    // Jeżeli łańcuch jest krótszy niż k, nie mamy co wycinać
    if (static_cast<int>(dna.size()) < k) {
        return spectrum;  // pusty wektor
    }

    // Wycinamy kolejne fragmenty o długości k
    for (int i = 0; i + k <= static_cast<int>(dna.size()); ++i) {
        spectrum.push_back(dna.substr(i, k));
    }
    return spectrum;
}