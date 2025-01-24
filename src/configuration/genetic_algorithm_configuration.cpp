//
// Created by konrad_guest on 23/01/2025.
//
#include "configuration/genetic_algorithm_configuration.h"
#include "metaheuristics/stopping_criteria.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <mutex>
#include <shared_mutex>

/*
    Konstruktor prywatny: ustawiamy domyślne wartości.
*/
GAConfig::GAConfig()
{
    std::cout << "[GAConfig] Creating new GAConfig instance" << std::endl;
    // All initialization is now done through in-class initializers
}

void GAConfig::resetToDefaults()
{
    std::unique_lock<std::shared_mutex> lock(configMutex);
    
    // Reset all values to their in-class initializer defaults
    m_maxGenerations = 100;
    populationSize = 100;
    mutationRate = 0.15;
    replacementRatio = 0.7;
    crossoverProbability = 1.0;
    
    selectionMethod = "tournament";
    crossoverType = "order";
    mutationMethod = "point";
    replacementMethod = "partial";
    stoppingMethod = "maxGenerations";
    
    noImprovementGenerations = 30;
    tournamentSize = 3;
    timeLimitSeconds = 60;
    
    adaptiveParams = {0.7, 20, 5, 0.1};
    
    alpha = 0.7;
    beta = 0.3;
    
    k = 8;
    deltaK = 2;
    lNeg = 25;
    lPoz = 25;
    repAllowed = false;
    probablePositive = 0;
    
    m_globalBestFit = -std::numeric_limits<double>::infinity();
    
    std::cout << "[GAConfig] Reset to default values" << std::endl;
}

// Setter z klamrowaniem wartości replacementRatio
void GAConfig::setReplacementRatio(double ratio) {
    std::unique_lock<std::shared_mutex> lock(configMutex);
    if (ratio >= 0.0 && ratio <= 1.0) {
        replacementRatio = ratio;
    } else {
        std::cerr << "[GAConfig] Invalid replacement ratio: " << ratio << ", using default of 0.7" << std::endl;
        replacementRatio = 0.7;
    }
}


// ---- IMPLEMENTACJA GET/SET GlobalBestFitness ----

double GAConfig::getGlobalBestFitness() const {
    std::shared_lock<std::shared_mutex> lock(configMutex);
    return m_globalBestFit;
}



void GAConfig::setGlobalBestFitness(double fitness) {
    std::unique_lock<std::shared_mutex> lock(configMutex);
    if (fitness > m_globalBestFit) {
        m_globalBestFit = fitness;
    }
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
void GAConfig::loadFromFile(const std::string& filePath)
{
    std::unique_lock<std::shared_mutex> lock(configMutex);
    
    if (isInitialized.load() && filePath == lastLoadedConfig) {
        std::cout << "[GAConfig] Configuration already loaded from " << filePath << std::endl;
        return;
    }
    
    std::cout << "[GAConfig] Loading configuration from " << filePath << std::endl;
    
    std::ifstream in(filePath);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open config file: " + filePath);
    }
    
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == ';') continue;
        
        size_t eqPos = line.find('=');
        if (eqPos == std::string::npos) continue;
        
        std::string key = line.substr(0, eqPos);
        std::string value = line.substr(eqPos + 1);
        
        // Trim whitespace
        key.erase(0, key.find_first_not_of(" \t"));
        key.erase(key.find_last_not_of(" \t") + 1);
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        
        try {
            if (key == "maxGenerations") {
                setMaxGenerations(std::stoi(value));
            }
            else if (key == "populationSize") {
                setPopulationSize(std::stoi(value));
            }
            else if (key == "mutationRate") {
                mutationRate = std::stod(value);
            }
            else if (key == "replacementRatio") {
                setReplacementRatio(std::stod(value));
            }
            else if (key == "crossoverProbability") {
                crossoverProbability = std::stod(value);
            }
            else if (key == "selectionMethod") {
                selectionMethod = value;
            }
            else if (key == "crossoverType") {
                crossoverType = value;
            }
            else if (key == "mutationMethod") {
                mutationMethod = value;
            }
            else if (key == "replacementMethod") {
                replacementMethod = value;
            }
            else if (key == "stoppingMethod") {
                stoppingMethod = value;
            }
            else if (key == "noImprovementGenerations") {
                noImprovementGenerations = std::stoi(value);
            }
            else if (key == "tournamentSize") {
                tournamentSize = std::stoi(value);
            }
            else if (key == "timeLimitSeconds") {
                timeLimitSeconds = std::stoi(value);
            }
            else if (key == "fitnessType") {
                fitnessType = value;
            }
            else if (key == "alpha") {
                alpha = std::stod(value);
            }
            else if (key == "beta") {
                beta = std::stod(value);
            }
            else if (key == "k") {
                k = std::stoi(value);
            }
            else if (key == "deltaK") {
                deltaK = std::stoi(value);
            }
            else if (key == "lNeg") {
                lNeg = std::stoi(value);
            }
            else if (key == "lPoz") {
                lPoz = std::stoi(value);
            }
            else if (key == "repAllowed") {
                repAllowed = (value == "true" || value == "1");
            }
            else if (key == "probablePositive") {
                probablePositive = std::stoi(value);
            }
            else if (key == "adaptive.inertia") {
                adaptiveParams.inertia = std::stod(value);
            }
            else if (key == "adaptive.adaptationInterval") {
                adaptiveParams.adaptationInterval = std::stoi(value);
            }
            else if (key == "adaptive.minTrials") {
                adaptiveParams.minTrials = std::stoi(value);
            }
            else if (key == "adaptive.minProb") {
                adaptiveParams.minProb = std::stod(value);
            }
            else {
                std::cerr << "[GAConfig] Unknown configuration key: " << key << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "[GAConfig] Error parsing " << key << ": " << e.what() << std::endl;
            continue;
        }
    }
    
    lastLoadedConfig = filePath;
    isInitialized.store(true);
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

// Add TimeLimitStopping class definition
class TimeLimitStopping : public IStopping {
public:
    explicit TimeLimitStopping(int sec) : limitSeconds(sec) {
        start = std::chrono::steady_clock::now();
    }
    
    bool stop(const std::vector<std::shared_ptr<std::vector<int>>>& population,
              int generation,
              const DNAInstance& instance,
              std::shared_ptr<IFitness> fitness,
              std::shared_ptr<IRepresentation> representation) override {
        auto now = std::chrono::steady_clock::now();
        auto elapsedSec = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
        return (elapsedSec >= limitSeconds);
    }
private:
    int limitSeconds;
    std::chrono::steady_clock::time_point start;
};

std::shared_ptr<IStopping> GAConfig::getStopping() const
{
    std::shared_lock<std::shared_mutex> lock(configMutex);
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
        return std::make_shared<TimeLimitStopping>(timeLimitSeconds);
    }
    // Default:
    std::cout << "[GAConfig] Using default stopping criteria" << std::endl;
    auto stopping = std::make_shared<MaxGenerationsStopping>(*const_cast<GAConfig*>(this));
    std::cout << "[GAConfig] Created default MaxGenerationsStopping using config value" << std::endl;
    return stopping;
}

void GAConfig::setParameters(const ParameterSet& ps) {
    std::unique_lock<std::shared_mutex> lock(configMutex);
    
    // Set parameters from the parameter set
    if (ps.params.count("maxGenerations")) {
        setMaxGenerations(std::stoi(ps.params.at("maxGenerations")));
    }
    if (ps.params.count("populationSize")) {
        setPopulationSize(std::stoi(ps.params.at("populationSize")));
    }
    if (ps.params.count("mutationRate")) {
        setMutationRate(std::stod(ps.params.at("mutationRate")));
    }
    if (ps.params.count("crossoverProbability")) {
        setCrossoverProbability(std::stod(ps.params.at("crossoverProbability")));
    }
    if (ps.params.count("replacementRatio")) {
        setReplacementRatio(std::stod(ps.params.at("replacementRatio")));
    }
    
    // Validate after setting all parameters
    validate();
}

bool GAConfig::validate() const {
    std::shared_lock<std::shared_mutex> lock(configMutex);
    
    bool isValid = true;
    
    if (m_maxGenerations <= 0) {
        std::cerr << "[GAConfig] Invalid maxGenerations: " << m_maxGenerations << std::endl;
        isValid = false;
    }
    
    if (populationSize <= 1) {
        std::cerr << "[GAConfig] Invalid population size: " << populationSize << std::endl;
        isValid = false;
    }
    
    if (mutationRate < 0.0 || mutationRate > 1.0) {
        std::cerr << "[GAConfig] Invalid mutation rate: " << mutationRate << std::endl;
        isValid = false;
    }
    
    if (crossoverProbability < 0.0 || crossoverProbability > 1.0) {
        std::cerr << "[GAConfig] Invalid crossover probability: " << crossoverProbability << std::endl;
        isValid = false;
    }
    
    if (replacementRatio < 0.0 || replacementRatio > 1.0) {
        std::cerr << "[GAConfig] Invalid replacement ratio: " << replacementRatio << std::endl;
        isValid = false;
    }
    
    return isValid;
}

