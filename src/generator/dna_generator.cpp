//
// Created by konrad_guest on 28/12/2024.
//

#include "generator/dna_generator.h"



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
    // by nie było 2-3 identycznych znaków pod rząd, itp. (wedle uznania).
    // Poniżej najprostsze rozwiązanie bez rozbudowanych reguł:
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
        // Jeżeli łańcuch jest krótszy niż k, nic nie tworzymy
        return spectrum;
    }

    auto &rng = RandomGenerator::getInstance().get();
    std::uniform_int_distribution<int> distDelta(0, deltaK);
    std::uniform_int_distribution<int> distSign(0, 1); // 0 => minus, 1 => plus

    // Pierwszy oligo zawsze ma długość k
    int startPos = 0;
    int currentOligoLength = k;

    // Ile w sumie powstanie oligonukleotydów = dna.size() - k + 1, ale z uwzględnieniem,
    // że "okno" przesuwamy o 1 za każdym razem.
    // Natomiast ostatnich (k+2) pozycji ma zawsze długość k.
    // Przykładowo – jeśli DNA ma 100 bp i k=8, to maksymalny indeks startowy to 92.
    // Ostatnie k+2 = 10 pozycji, więc od 90 do 99 włącznie.
    // Dlatego do momentu ( dna.size() - (k+2) ) wylosujemy ewentualnie inne długości.

    int endFixedZoneStart = dna.size() - (k + 2);
    if(endFixedZoneStart < 0) {
        // Gdyby k+2 > dna.size(), wówczas nie ma "zmiennej" części, wszystko jest stałe.
        endFixedZoneStart = 0;
    }

    while (startPos + k <= (int)dna.size()) {
        // Sprawdzamy, czy jesteśmy już w strefie ostatnich (k+2) oligonów
        if (startPos <= endFixedZoneStart) {
            // Dla pierwszego (startPos == 0) zachowujemy k
            if (startPos == 0) {
                currentOligoLength = k;
            } else {
                // Losujemy wartość z zakresu [0..deltaK]
                int d = distDelta(rng);
                if (d > 0) {
                    // Losujemy znak
                    int sign = distSign(rng);
                    if (sign == 0) {
                        currentOligoLength = std::max(k - d, 1); 
                    } else {
                        currentOligoLength = std::min(k + d, (int)dna.size() - startPos);
                    }
                } else {
                    // D == 0 => nic nie zmieniamy
                    currentOligoLength = k;
                }
            }
        } else {
            // W strefie końcowej zawsze bierzemy oligo o długości k,
            // o ile wystarczy znaków do wzięcia
            if (startPos + k > (int)dna.size()) {
                break;
            }
            currentOligoLength = k;
        }

        // Zabezpieczenie, żeby nie wyjść za koniec łańcucha
        if (startPos + currentOligoLength > (int)dna.size()) {
            break;
        }

        // Wycinamy fragment
        spectrum.push_back(dna.substr(startPos, currentOligoLength));

        // Przesuwamy okno o 1
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

    // Zapisujemy: n, k, deltaK, lNeg, lPoz, repAllowed, probablePositive
    out << instance.getN() << " " 
        << instance.getK() << " "
        << instance.getDeltaK() << " "
        << instance.getLNeg() << " "
        << instance.getLPoz() << " "
        << (instance.isRepAllowed() ? 1 : 0) << " "
        << instance.getProbablePositive() << "\n";

    // Zapisujemy DNA
    out << instance.getDNA() << "\n";

    // Zapisujemy rozmiar spektrum i poszczególne k-mery
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
    // Walidacja, np. min = 300, max = 700
    if(n < 300 || n > 700) {
        // Można rzucić wyjątek, lub ustawić na domyślną wartość
        n = 400;
    }
    instance.setN(n);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setK(int k)
{
    // Walidacja, np. min = 7, max = 10
    if(k < 7 || k > 10) {
        k = 8;
    }
    instance.setK(k);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setDeltaK(int dk)
{
    // Walidacja, np. min = 0, max = 2
    if(dk < 0 || dk > 2) {
        dk = 2;
    }
    instance.setDeltaK(dk);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLNeg(int ln)
{
    // l_neg może być 0 lub minimum 10,
    // nie więcej niż ~15% wszystkich oligonukleotydów (to sprawdza się dopiero po wytworzeniu?)
    // Na potrzeby przykładu – wstępna walidacja:
    if(ln != 0 && ln < 10) {
        ln = 0;  // lub throw
    }
    instance.setLNeg(ln);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setLPoz(int lp) {
    // Sprawdzamy czy liczba błędów pozytywnych jest rozsądna
    if(lp != 0) {
        // Minimum 10 błędów jeśli nie zero
        if(lp < 10) lp = 10;
        
        // Ograniczamy maksymalną liczbę błędów do ~50% obecnego spektrum
        // aby uniknąć zbyt dużego rozrostu
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
    // Parametr ignorowany, jeśli l_neg > 0
    if(instance.getLNeg() > 0) {
        rep = true;  // lub false – zależnie od interpretacji; w treści "nieistotne"
    }
    instance.setRepAllowed(rep);
    return *this;
}

DNAInstanceBuilder& DNAInstanceBuilder::setProbablePositive(int val)
{
    // 0 lub 1
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

    // Tutaj można dodatkowo sprawdzić czy lNeg <= 15% wielkości spektrum, itp.
    // lub ewentualnie skorygować te wartości.

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

    // Nie można usunąć więcej niż istnieje
    int toRemove = std::min((int)spectrum.size(), lNeg);

    auto &rng = RandomGenerator::getInstance().get();
    std::shuffle(spectrum.begin(), spectrum.end(), rng);

    spectrum.erase(spectrum.begin(), spectrum.begin() + toRemove);
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
    
    // Helper do generowania losowego k-meru
    auto generateRandomKmer = [&](int length) {
        std::string kmer;
        kmer.reserve(length);
        for(int i=0; i<length; i++) {
            kmer.push_back(nucleotides[nucDist(rng)]);
        }
        return kmer;
    };
    
    // Dodajemy błędy pozytywne
    int errorsAdded = 0;
    int maxAttempts = lPoz * 10; // Limit prób, aby uniknąć nieskończonej pętli
    int attempts = 0;
    
    while(errorsAdded < lPoz && attempts < maxAttempts) {
        attempts++;
        
        if(probablePositive == 0) {
            // Generowanie losowej sekwencji
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
            // Modyfikacja istniejącego oligonukleotydu
            std::uniform_int_distribution<int> indexDist(0, spectrum.size() - 1);
            int idx = indexDist(rng);
            std::string original = spectrum[idx];
            
            // Modyfikujemy kopię zmieniając jeden losowy nukleotyd
            std::string modified = original;
            int posToChange = (distSign(rng) == 0) ? 
                modified.length() - 1 :  // ostatni
                modified.length() / 2;   // środkowy
                
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
        std::cerr << "Warning: Could only add " << errorsAdded 
                  << " positive errors out of " << lPoz << " requested\n";
    }
}