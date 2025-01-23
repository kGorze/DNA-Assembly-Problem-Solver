//
// Created by konrad_guest on 23/01/2025.
//
#include "configuration/genetic_algorithm_configuration.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <chrono>

/*
    Konstruktor prywatny: ustawiamy domyślne wartości.
*/
GAConfig::GAConfig()
{
    // Domyślne wartości
    populationSize       = 100;
    mutationRate         = 0.15;
    replacementRatio     = 0.7;
    maxGenerations       = 200;
    crossoverProbability = 1.0; // zawsze krzyżujemy

    selectionMethod      = "tournament";
    crossoverType        = "order";
    mutationMethod       = "point";
    replacementMethod    = "partial";
    stoppingMethod       = "maxGenerations";

    noImprovementGenerations = 30;
    tournamentSize           = 3;
    timeLimitSeconds         = 60; // np. 60s

    adaptiveParams.inertia            = 0.7;
    adaptiveParams.adaptationInterval = 20;
    adaptiveParams.minTrials          = 5;
    adaptiveParams.minProb            = 0.1;

    fitnessType = "optimized_graph";
    alpha = 0.7;
    beta  = 0.3;

    // Najlepszy fitness domyślnie -∞

    m_globalBestFit = -std::numeric_limits<double>::infinity();
}

GAConfig& GAConfig::getInstance()
{
    static GAConfig instance;
    return instance;
}

// Setter z klamrowaniem wartości replacementRatio
void GAConfig::setReplacementRatio(double ratio) {
    if (ratio < 0.0) ratio = 0.0;
    if (ratio > 1.0) ratio = 1.0;
    replacementRatio = ratio;
}


// ---- IMPLEMENTACJA GET/SET GlobalBestFitness ----

double GAConfig::getGlobalBestFitness() const {

    return m_globalBestFit;

}



void GAConfig::setGlobalBestFitness(double fitness) {

    m_globalBestFit = fitness;

}


/*
    Prosty parser pliku konfiguracyjnego w formacie:
    -----------------------------------
    # Komentarz
    populationSize=200
    mutationRate=0.25
    selectionMethod=tournament
    ...
    # Komentarz
    adaptive.inertia=0.6
    adaptive.adaptationInterval=15
    ...
    -----------------------------------

    - Linie zaczynające się od '#' lub ';' traktujemy jako komentarz.
    - Puste linie ignorujemy.
    - Klucz i wartość rozdzielone '='.
    - Dla prostoty: "adaptive.xyz" -> interpretujemy i przypisujemy do struct AdaptiveCrossoverParams.

*/
bool GAConfig::loadFromFile(const std::string& filePath)
{
    std::vector<std::string> possiblePaths = {
        filePath,                // Oryginalnie podana ścieżka
        "config.cfg",
        "./config.cfg",
        "../config.cfg",
        "../build/config.cfg",
        "./build/config.cfg"// Dodana ścieżka do folderu build
    };

    for (const auto& path : possiblePaths) {
        std::ifstream in(path);
        if (in.is_open()) {
            std::cout << "[GAConfig] Ładowanie pliku konfiguracyjnego z: " << path << std::endl;

            std::string line;
            while (std::getline(in, line)) {
                // Usuwamy spacje z przodu i końca
                auto trim = [&](std::string &s) {
                    // Trim z lewej strony
                    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch){
                        return !std::isspace(ch);
                    }));
                    // Trim z prawej strony
                    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch){
                        return !std::isspace(ch);
                    }).base(), s.end());
                };
                trim(line);

                // Pomijamy komentarze (#, ;) i puste linie
                if (line.empty() || line[0] == '#' || line[0] == ';')
                    continue;

                // Szukamy '='
                auto eqPos = line.find('=');
                if (eqPos == std::string::npos) {
                    // Brak '=' -> pomijamy
                    continue;
                }

                std::string key   = line.substr(0, eqPos);
                std::string value = line.substr(eqPos + 1);
                trim(key);
                trim(value);

                // Parsujemy klucz i przypisujemy wartości do odpowiednich pól
                try {
                    if (key == "populationSize") {
                        populationSize = std::stoi(value);
                    } else if (key == "mutationRate") {
                        mutationRate = std::stod(value);
                    } else if (key == "replacementRatio") {
                        setReplacementRatio(std::stod(value));
                    } else if (key == "maxGenerations") {
                        maxGenerations = std::stoi(value);
                    } else if (key == "crossoverProbability") {
                        crossoverProbability = std::stod(value);
                    } else if (key == "selectionMethod") {
                        selectionMethod = value;
                    } else if (key == "crossoverType") {
                        crossoverType = value;
                    } else if (key == "mutationMethod") {
                        mutationMethod = value;
                    } else if (key == "replacementMethod") {
                        replacementMethod = value;
                    } else if (key == "stoppingMethod") {
                        stoppingMethod = value;
                    } else if (key == "noImprovementGenerations") {
                        noImprovementGenerations = std::stoi(value);
                    } else if (key == "tournamentSize") {
                        tournamentSize = std::stoi(value);
                    } else if (key == "timeLimitSeconds") {
                        timeLimitSeconds = std::stoi(value);
                    } else if (key == "fitnessType") {
                        fitnessType = value;
                    } else if (key == "alpha") {
                        alpha = std::stod(value);
                    } else if (key == "beta") {
                        beta = std::stod(value);
                    
                    // Obsługa parametrów adaptacyjnego krzyżowania (przykład):
                    } else if (key == "adaptive.inertia") {
                        adaptiveParams.inertia = std::stod(value);
                    } else if (key == "adaptive.adaptationInterval") {
                        adaptiveParams.adaptationInterval = std::stoi(value);
                    } else if (key == "adaptive.minTrials") {
                        adaptiveParams.minTrials = std::stoi(value);
                    } else if (key == "adaptive.minProb") {
                        adaptiveParams.minProb = std::stod(value);
                    } else {
                        // Nieznany klucz – wypisz ostrzeżenie
                        std::cerr << "[GAConfig] Nieznany klucz konfiguracyjny: " << key << std::endl;
                    }
                } catch (const std::exception& e) {
                    std::cerr << "[GAConfig] Błąd parsowania wartości dla klucza '" << key << "': " << e.what() << std::endl;
                }
            }

            // Zamykamy plik po zakończeniu parsowania
            in.close();
            return true; // Sukces – plik został załadowany
        }
    }

    // Jeśli żaden plik nie został znaleziony i załadowany
    std::cerr << "[GAConfig] Nie znaleziono żadnego pliku konfiguracyjnego w dostępnych ścieżkach." << std::endl;
    return false;
}

// ---------------------------------------
// Metody zwracające obiekty metaheurystyk
// ---------------------------------------
std::shared_ptr<IRepresentation> GAConfig::getRepresentation() const
{
    // Na potrzeby SBH zwykle PermutationRepresentation,
    // można tu dodać logikę wyboru różnych representation (jeśli potrzeba).
    return std::make_shared<PermutationRepresentation>();
}

std::shared_ptr<ISelection> GAConfig::getSelection() const
{
    // W zależności od selectionMethod
    if (selectionMethod == "tournament") {
        return std::make_shared<TournamentSelection>(tournamentSize, m_cache);
    }
    // else if (selectionMethod == "roulette") { ... }

    // Domyślnie:
    return std::make_shared<TournamentSelection>(tournamentSize, m_cache);
}

std::shared_ptr<ICrossover> GAConfig::getCrossover(const std::string& type) const
{
    if (type == "order") {
        return std::make_shared<OrderCrossover>();
    } else if (type == "edge") {
        return std::make_shared<EdgeRecombination>();
    } else if (type == "pmx") {
        return std::make_shared<PMXCrossover>();
    } else if (type == "onepoint") {
        return std::make_shared<OnePointCrossover>();
    } else if (type == "adaptive") {
        auto ac = std::make_shared<AdaptiveCrossover>();
        ac->setParameters(adaptiveParams.inertia,
                          adaptiveParams.adaptationInterval,
                          adaptiveParams.minTrials,
                          adaptiveParams.minProb);
        return ac;
    }
    // Domyślnie:
    return std::make_shared<OrderCrossover>();
}

std::shared_ptr<IMutation> GAConfig::getMutation() const
{
    if (mutationMethod == "point") {
        return std::make_shared<PointMutation>(mutationRate);
    }
    // Możesz zdefiniować inne klasy:
    //  if (mutationMethod == "swap") { ... }
    //  if (mutationMethod == "scramble") { ... }
    // itp.

    // Domyślnie point:
    return std::make_shared<PointMutation>(mutationRate);
}

std::shared_ptr<IReplacement> GAConfig::getReplacement() const
{
    if (replacementMethod == "full") {
        return std::make_shared<FullReplacement>();
    } else if (replacementMethod == "partial") {
        return std::make_shared<PartialReplacement>(replacementRatio, m_cache);
    }
    // domyślnie partial
    return std::make_shared<PartialReplacement>(replacementRatio, m_cache);
}

std::shared_ptr<IFitness> GAConfig::getFitness() const
{
    if (fitnessType == "simple") {
        return std::make_shared<SimpleFitness>();
    } else if (fitnessType == "better") {
        return std::make_shared<BetterFitness>();
    } else if (fitnessType == "smithwaterman") {
        return std::make_shared<SmithWatermanFitness>();
    } else if (fitnessType == "optimized_graph") {
        // ewentualnie można wstrzykiwać alpha, beta, jeśli zdefiniujemy settery w tej klasie
        auto of = std::make_shared<OptimizedGraphBasedFitness>();
        return of;
    }
    // Domyślnie:
    return std::make_shared<BetterFitness>();
}

std::shared_ptr<IStopping> GAConfig::getStopping() const
{
    if (stoppingMethod == "maxGenerations") {
        return std::make_shared<MaxGenerationsStopping>(maxGenerations);
    } else if (stoppingMethod == "noImprovement") {
        return std::make_shared<NoImprovementStopping>(noImprovementGenerations);
    } else if (stoppingMethod == "timeLimit") {
        // Przykład prostej implementacji TimeLimitStopping
        class TimeLimitStopping : public IStopping {
        public:
            TimeLimitStopping(int sec) : limitSeconds(sec) {
                start = std::chrono::steady_clock::now();
            }
            bool stop(const std::vector<std::shared_ptr<std::vector<int>>> &population,
                      int generation,
                      const DNAInstance &instance,
                      std::shared_ptr<IFitness> fitness,
                      std::shared_ptr<IRepresentation> representation) override
            {
                auto now = std::chrono::steady_clock::now();
                auto elapsedSec = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
                return (elapsedSec >= limitSeconds);
            }
        private:
            int limitSeconds;
            std::chrono::steady_clock::time_point start;
        };

        return std::make_shared<TimeLimitStopping>(timeLimitSeconds);
    }
    // Domyślnie:
    return std::make_shared<MaxGenerationsStopping>(maxGenerations);
}

