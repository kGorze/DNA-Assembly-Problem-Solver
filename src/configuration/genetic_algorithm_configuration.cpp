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
#include <mutex>

/*
    Konstruktor prywatny: ustawiamy domyślne wartości.
*/
GAConfig::GAConfig()
{
    std::cout << "[GAConfig] Creating new GAConfig instance" << std::endl;
    // Don't reset to defaults automatically
    m_maxGenerations = -1;  // Invalid value to indicate it needs to be set
    isInitialized = false;
}

void GAConfig::resetToDefaults() 
{
    std::lock_guard<std::mutex> lock(configMutex);
    std::cout << "[GAConfig] Resetting to default values" << std::endl;
    
    // Only set defaults if not already initialized from config
    if (!isInitialized) {
        // Domyślne wartości
        populationSize       = 100;
        mutationRate         = 0.15;
        replacementRatio     = 0.7;
        m_maxGenerations     = 200;  // Using private member
        crossoverProbability = 1.0;

        selectionMethod      = "tournament";
        crossoverType        = "order";
        mutationMethod       = "point";
        replacementMethod    = "partial";
        stoppingMethod       = "maxGenerations";

        noImprovementGenerations = 30;
        tournamentSize           = 3;
        timeLimitSeconds         = 60;

        adaptiveParams.inertia            = 0.7;
        adaptiveParams.adaptationInterval = 20;
        adaptiveParams.minTrials          = 5;
        adaptiveParams.minProb            = 0.1;

        fitnessType = "optimized_graph";
        alpha = 0.7;
        beta  = 0.3;

        // DNA Generation parameters - invalid values to indicate they need to be set
        k = -1;
        deltaK = -1;
        lNeg = -1;
        lPoz = -1;
        repAllowed = false;
        probablePositive = -1;

        m_globalBestFit = -std::numeric_limits<double>::infinity();
        
        std::cout << "[GAConfig] Default values set:" << std::endl
                  << "  maxGenerations = " << m_maxGenerations << std::endl
                  << "  populationSize = " << populationSize << std::endl
                  << "  mutationRate = " << mutationRate << std::endl;
    } else {
        std::cout << "[GAConfig] Already initialized from config, not resetting" << std::endl;
    }
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
    std::lock_guard<std::mutex> lock(configMutex);
    
    std::cout << "[GAConfig] Current maxGenerations before loading: " << m_maxGenerations << std::endl;
    
    if (filePath == lastLoadedConfig && isInitialized) {
        std::cout << "[GAConfig] Config file " << filePath << " was already loaded, skipping." << std::endl;
        return true;
    }
    
    // Get absolute path of executable
    char exePath[1024];
    #ifdef _WIN32
    GetModuleFileName(NULL, exePath, sizeof(exePath));
    #else
    if (readlink("/proc/self/exe", exePath, sizeof(exePath)) == -1) {
        LOG_ERROR("Failed to get executable path");
        return false;
    }
    #endif
    std::string exeDir = std::string(exePath);
    exeDir = exeDir.substr(0, exeDir.find_last_of("/\\"));
    
    std::vector<std::string> possiblePaths = {
        filePath,
        "config.cfg",
        "./config.cfg",
        "../config.cfg",
        (exeDir + "/config.cfg"),
        (exeDir + "/../config.cfg"),
        (exeDir + "/build/config.cfg"),
        "../build/config.cfg",
        "./build/config.cfg"
    };

    for (const auto& path : possiblePaths) {
        std::ifstream in(path);
        if (in.is_open()) {
            std::cout << "[GAConfig] Loading configuration from: " << path << std::endl;

            std::string line;
            bool anyValueSet = false;
            
            while (std::getline(in, line)) {
                auto trim = [&](std::string &s) {
                    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch){
                        return !std::isspace(ch);
                    }));
                    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch){
                        return !std::isspace(ch);
                    }).base(), s.end());
                };
                trim(line);

                if (line.empty() || line[0] == '#' || line[0] == ';') {
                    continue;
                }

                auto eqPos = line.find('=');
                if (eqPos == std::string::npos) {
                    continue;
                }

                std::string key = line.substr(0, eqPos);
                std::string value = line.substr(eqPos + 1);
                trim(key);
                trim(value);

                std::cout << "[GAConfig] Processing: " << key << " = " << value << std::endl;

                try {
                    if (key == "maxGenerations") {
                        int oldValue = m_maxGenerations;
                        m_maxGenerations = std::stoi(value);
                        std::cout << "[GAConfig] Updated maxGenerations from " << oldValue << " to " << m_maxGenerations << std::endl;
                        anyValueSet = true;
                    } else if (key == "populationSize") {
                        populationSize = std::stoi(value);
                        std::cout << "[GAConfig] Set populationSize = " << populationSize << std::endl;
                    } else if (key == "mutationRate") {
                        mutationRate = std::stod(value);
                        std::cout << "[GAConfig] Set mutationRate = " << mutationRate << std::endl;
                    } else if (key == "replacementRatio") {
                        setReplacementRatio(std::stod(value));
                        std::cout << "[GAConfig] Set replacementRatio = " << replacementRatio << std::endl;
                    } else if (key == "crossoverProbability") {
                        crossoverProbability = std::stod(value);
                        std::cout << "[GAConfig] Set crossoverProbability = " << crossoverProbability << std::endl;
                    } else if (key == "selectionMethod") {
                        selectionMethod = value;
                        std::cout << "[GAConfig] Set selectionMethod = " << selectionMethod << std::endl;
                    } else if (key == "crossoverType") {
                        crossoverType = value;
                        std::cout << "[GAConfig] Set crossoverType = " << crossoverType << std::endl;
                    } else if (key == "mutationMethod") {
                        mutationMethod = value;
                        std::cout << "[GAConfig] Set mutationMethod = " << mutationMethod << std::endl;
                    } else if (key == "replacementMethod") {
                        replacementMethod = value;
                        std::cout << "[GAConfig] Set replacementMethod = " << replacementMethod << std::endl;
                    } else if (key == "stoppingMethod") {
                        stoppingMethod = value;
                        std::cout << "[GAConfig] Set stoppingMethod = " << stoppingMethod << std::endl;
                    } else if (key == "noImprovementGenerations") {
                        noImprovementGenerations = std::stoi(value);
                        std::cout << "[GAConfig] Set noImprovementGenerations = " << noImprovementGenerations << std::endl;
                    } else if (key == "tournamentSize") {
                        tournamentSize = std::stoi(value);
                        std::cout << "[GAConfig] Set tournamentSize = " << tournamentSize << std::endl;
                    } else if (key == "timeLimitSeconds") {
                        timeLimitSeconds = std::stoi(value);
                        std::cout << "[GAConfig] Set timeLimitSeconds = " << timeLimitSeconds << std::endl;
                    } else if (key == "fitnessType") {
                        fitnessType = value;
                        std::cout << "[GAConfig] Set fitnessType = " << fitnessType << std::endl;
                    } else if (key == "alpha") {
                        alpha = std::stod(value);
                        std::cout << "[GAConfig] Set alpha = " << alpha << std::endl;
                    } else if (key == "beta") {
                        beta = std::stod(value);
                        std::cout << "[GAConfig] Set beta = " << beta << std::endl;
                    
                    // Instance parameters
                    } else if (key == "k") {
                        k = std::stoi(value);
                        std::cout << "[GAConfig] Set k = " << k << std::endl;
                    } else if (key == "deltaK") {
                        deltaK = std::stoi(value);
                        std::cout << "[GAConfig] Set deltaK = " << deltaK << std::endl;
                    } else if (key == "lNeg") {
                        lNeg = std::stoi(value);
                        std::cout << "[GAConfig] Set lNeg = " << lNeg << std::endl;
                    } else if (key == "lPoz") {
                        lPoz = std::stoi(value);
                        std::cout << "[GAConfig] Set lPoz = " << lPoz << std::endl;
                    } else if (key == "repAllowed") {
                        repAllowed = (value == "true" || value == "1");
                        std::cout << "[GAConfig] Set repAllowed = " << repAllowed << std::endl;
                    } else if (key == "probablePositive") {
                        probablePositive = std::stoi(value);
                        std::cout << "[GAConfig] Set probablePositive = " << probablePositive << std::endl;
                    
                    // Obsługa parametrów adaptacyjnego krzyżowania (przykład):
                    } else if (key == "adaptive.inertia") {
                        adaptiveParams.inertia = std::stod(value);
                        std::cout << "[GAConfig] Set adaptive.inertia = " << adaptiveParams.inertia << std::endl;
                    } else if (key == "adaptive.adaptationInterval") {
                        adaptiveParams.adaptationInterval = std::stoi(value);
                        std::cout << "[GAConfig] Set adaptive.adaptationInterval = " << adaptiveParams.adaptationInterval << std::endl;
                    } else if (key == "adaptive.minTrials") {
                        adaptiveParams.minTrials = std::stoi(value);
                        std::cout << "[GAConfig] Set adaptive.minTrials = " << adaptiveParams.minTrials << std::endl;
                    } else if (key == "adaptive.minProb") {
                        adaptiveParams.minProb = std::stod(value);
                        std::cout << "[GAConfig] Set adaptive.minProb = " << adaptiveParams.minProb << std::endl;
                    } else {
                        // Nieznany klucz – wypisz ostrzeżenie
                        std::cerr << "[GAConfig] Unknown configuration key: " << key << std::endl;
                    }
                } catch (const std::exception& e) {
                    std::cerr << "[GAConfig] Error parsing value for key '" << key << "': " << e.what() << std::endl;
                }
            }

            if (anyValueSet) {
                lastLoadedConfig = path;
                isInitialized = true;
            }

            std::cout << "\n[GAConfig] Final configuration values:" << std::endl
                      << "  maxGenerations = " << m_maxGenerations << std::endl
                      << "  populationSize = " << populationSize << std::endl
                      << "  mutationRate = " << mutationRate << std::endl;

            in.close();
            return true;
        }
    }

    std::cerr << "[GAConfig] No configuration file found in available paths." << std::endl;
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
        return std::make_shared<TournamentSelection>(*const_cast<GAConfig*>(this), m_cache);
    }
    // else if (selectionMethod == "roulette") { ... }

    // Domyślnie:
    return std::make_shared<TournamentSelection>(*const_cast<GAConfig*>(this), m_cache);
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
    std::lock_guard<std::mutex> lock(configMutex);
    // Log the current maxGenerations value for debugging
    std::cout << "[GAConfig] Creating stopping criteria with maxGenerations=" << m_maxGenerations << std::endl;
    
    if (stoppingMethod == "maxGenerations") {
        // Pass this instance to the constructor
        auto stopping = std::make_shared<MaxGenerationsStopping>(*const_cast<GAConfig*>(this));
        std::cout << "[GAConfig] Created MaxGenerationsStopping using config value" << std::endl;
        return stopping;
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
    // Default:
    std::cout << "[GAConfig] Using default stopping criteria" << std::endl;
    auto stopping = std::make_shared<MaxGenerationsStopping>(*const_cast<GAConfig*>(this));
    std::cout << "[GAConfig] Created default MaxGenerationsStopping using config value" << std::endl;
    return stopping;
}

void GAConfig::setParameters(const ParameterSet &ps) {
    std::lock_guard<std::mutex> lock(configMutex);
    
    if (ps.params.count("populationSize")) {
        populationSize = std::stoi(ps.params.at("populationSize"));
    }
    if (ps.params.count("mutationRate")) {
        mutationRate = std::stod(ps.params.at("mutationRate"));
    }
    if (ps.params.count("replacementRatio")) {
        replacementRatio = std::stod(ps.params.at("replacementRatio"));
    }
    if (ps.params.count("tournamentSize")) {
        tournamentSize = std::stoi(ps.params.at("tournamentSize"));
    }
    if (ps.params.count("crossoverType")) {
        crossoverType = ps.params.at("crossoverType");
    }
    if (ps.params.count("selectionMethod")) {
        selectionMethod = ps.params.at("selectionMethod");
    }
    if (ps.params.count("adaptive.inertia")) {
        adaptiveParams.inertia = std::stod(ps.params.at("adaptive.inertia"));
    }
    if (ps.params.count("adaptive.adaptationInterval")) {
        adaptiveParams.adaptationInterval = std::stoi(ps.params.at("adaptive.adaptationInterval"));
    }
    if (ps.params.count("adaptive.minTrials")) {
        adaptiveParams.minTrials = std::stoi(ps.params.at("adaptive.minTrials"));
    }
    if (ps.params.count("adaptive.minProb")) {
        adaptiveParams.minProb = std::stod(ps.params.at("adaptive.minProb"));
    }
    if (ps.params.count("k")) {
        k = std::stoi(ps.params.at("k"));
    }
    if (ps.params.count("deltaK")) {
        deltaK = std::stoi(ps.params.at("deltaK"));
    }
    if (ps.params.count("lNeg")) {
        lNeg = std::stoi(ps.params.at("lNeg"));
    }
    if (ps.params.count("lPoz")) {
        lPoz = std::stoi(ps.params.at("lPoz"));
    }
    if (ps.params.count("repAllowed")) {
        repAllowed = ps.params.at("repAllowed") == "true";
    }
    if (ps.params.count("probablePositive")) {
        probablePositive = std::stoi(ps.params.at("probablePositive"));
    }
}

bool GAConfig::validate() const {
    if (k <= 0) {
        std::cerr << "[ERROR] K must be positive" << std::endl;
        return false;
    }
    
    if (populationSize <= 0) {
        std::cerr << "[ERROR] Population size must be positive" << std::endl;
        return false;
    }
    
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        std::cerr << "[ERROR] Mutation rate must be between 0 and 1" << std::endl;
        return false;
    }
    
    if (replacementRatio < 0.0 || replacementRatio > 1.0) {
        std::cerr << "[ERROR] Replacement ratio must be between 0 and 1" << std::endl;
        return false;
    }
    
    return true;
}

