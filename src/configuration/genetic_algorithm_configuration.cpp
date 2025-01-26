//
// Created by konrad_guest on 23/01/2025.
//
#include "../../include/configuration/genetic_algorithm_configuration.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <chrono>
#include <mutex>
#include <shared_mutex>
#include <limits>
#include "../../include/metaheuristics/selection_impl.h"
#include "../../include/metaheuristics/crossover_impl.h"
#include "../../include/metaheuristics/mutation_impl.h"
#include "../../include/metaheuristics/replacement_impl.h"
#include "../../include/metaheuristics/fitness_impl.h"
#include "../../include/metaheuristics/stopping_criteria_impl.h"
#include "../../include/interfaces/i_population_cache.h"
#include "../../include/metaheuristics/population_cache_impl.h"
#include "../../include/metaheuristics/representation_impl.h"
#include "../../include/metaheuristics/adaptive_crossover.h"
#include <json/json.h>

/*
    Konstruktor prywatny: ustawiamy domyślne wartości.
*/
// GAConfig::GAConfig() { ... }

void GAConfig::resetToDefaults()
{
    std::unique_lock<std::shared_mutex> lock(configMutex);
    
    m_maxGenerations = 1000;
    m_populationSize = 100;
    m_mutationRate = 0.1;
    replacementRatio = 0.8;
    crossoverProbability = 0.8;
    targetFitness = 1.0;
    
    selectionMethod = "tournament";
    crossoverType = "order";
    mutationMethod = "point";
    replacementMethod = "partial";
    stoppingMethod = "maxGenerations";
    fitnessType = "optimized_graph";
    
    noImprovementGenerations = 30;
    tournamentSize = 5;
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

bool GAConfig::loadFromFile(const std::string& filePath)
{
    // First check if we've already loaded this config
    if (isInitialized.load() && filePath == lastLoadedConfig) {
        std::cout << "[GAConfig] Configuration already loaded from " << filePath << std::endl;
        return true;
    }
    
    std::cout << "[GAConfig] Loading configuration from " << filePath << std::endl;
    
    // Open and read the file without holding the lock
    std::ifstream in(filePath);
    if (!in.is_open()) {
        std::string error = "Failed to open config file: " + filePath;
        throw std::runtime_error(error);
    }
    
    // Create temporary storage for the new values
    int maxGen = 1000;
    int popSize = 100;
    double mutRate = 0.1;
    double replRatio = 0.8;
    double crossProb = 0.8;
    std::string selMethod = "tournament";
    std::string crossType = "order";
    std::string mutMethod = "point";
    std::string replMethod = "partial";
    std::string stopMethod = "maxGenerations";
    std::string fitType = "optimized_graph";
    int noImprovGen = 30;
    int tournSize = 5;
    int timeLimit = 60;
    AdaptiveCrossoverParams adaptParams = {0.7, 20, 5, 0.1};
    double a = 0.7, b = 0.3;
    int k_val = 8, dk = 2, ln = 25, lp = 25;
    bool rep = false;
    int probPos = 0;
    double targetFitness = 1.0;
    
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
            if (key == "maxGenerations") maxGen = std::stoi(value);
            else if (key == "populationSize") popSize = std::stoi(value);
            else if (key == "mutationRate") mutRate = std::stod(value);
            else if (key == "replacementRatio") replRatio = std::stod(value);
            else if (key == "crossoverProbability") crossProb = std::stod(value);
            else if (key == "selectionMethod") selMethod = value;
            else if (key == "crossoverType") crossType = value;
            else if (key == "mutationMethod") mutMethod = value;
            else if (key == "replacementMethod") replMethod = value;
            else if (key == "stoppingMethod") stopMethod = value;
            else if (key == "fitnessType") fitType = value;
            else if (key == "noImprovementGenerations") noImprovGen = std::stoi(value);
            else if (key == "tournamentSize") tournSize = std::stoi(value);
            else if (key == "timeLimitSeconds") timeLimit = std::stoi(value);
            else if (key == "alpha") a = std::stod(value);
            else if (key == "beta") b = std::stod(value);
            else if (key == "k") k_val = std::stoi(value);
            else if (key == "deltaK") dk = std::stoi(value);
            else if (key == "lNeg") ln = std::stoi(value);
            else if (key == "lPoz") lp = std::stoi(value);
            else if (key == "repAllowed") rep = (value == "true" || value == "1");
            else if (key == "probablePositive") probPos = std::stoi(value);
            else if (key == "adaptive.inertia") adaptParams.inertia = std::stod(value);
            else if (key == "adaptive.adaptationInterval") adaptParams.adaptationInterval = std::stoi(value);
            else if (key == "adaptive.minTrials") adaptParams.minTrials = std::stoi(value);
            else if (key == "adaptive.minProb") adaptParams.minProb = std::stod(value);
            else if (key == "targetFitness") targetFitness = std::stod(value);
        } catch (const std::exception& e) {
            std::cerr << "[GAConfig] Error parsing " << key << ": " << e.what() << std::endl;
        }
    }
    
    // Now acquire the lock and update all values atomically
    std::unique_lock<std::shared_mutex> lock(configMutex);
    
    m_maxGenerations = maxGen;
    m_populationSize = popSize;
    m_mutationRate = mutRate;
    replacementRatio = replRatio;
    crossoverProbability = crossProb;
    selectionMethod = selMethod;
    crossoverType = crossType;
    mutationMethod = mutMethod;
    replacementMethod = replMethod;
    stoppingMethod = stopMethod;
    fitnessType = fitType;
    noImprovementGenerations = noImprovGen;
    tournamentSize = tournSize;
    timeLimitSeconds = timeLimit;
    adaptiveParams = adaptParams;
    alpha = a;
    beta = b;
    k = k_val;
    deltaK = dk;
    lNeg = ln;
    lPoz = lp;
    repAllowed = rep;
    probablePositive = probPos;
    targetFitness = targetFitness;
    
    lastLoadedConfig = filePath;
    isInitialized.store(true);
    
    std::cout << "[GAConfig] Configuration loaded successfully" << std::endl;
    return true;
}

// ---------------------------------------
// Metody zwracające obiekty metaheurystyk
// ---------------------------------------
std::shared_ptr<IRepresentation> GAConfig::getRepresentation() const
{
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
        return std::make_shared<OnePointCrossover>(crossoverProbability);
    } else if (type == "adaptive") {
        auto ac = std::make_shared<AdaptiveCrossover>();
        ac->setParameters(adaptiveParams.inertia,
                          adaptiveParams.adaptationInterval,
                          adaptiveParams.minTrials,
                          adaptiveParams.minProb);
        return ac;
    }
    // Default:
    return std::make_shared<OrderCrossover>();
}

std::shared_ptr<IMutation> GAConfig::getMutation() const
{
    if (mutationMethod == "point") {
        return std::make_shared<PointMutation>(m_mutationRate);
    }
    // Możesz zdefiniować inne klasy:
    //  if (mutationMethod == "swap") { ... }
    //  if (mutationMethod == "scramble") { ... }
    // itp.

    // Domyślnie point:
    return std::make_shared<PointMutation>(m_mutationRate);
}

std::shared_ptr<IReplacement> GAConfig::getReplacement() const
{
    if (replacementMethod == "full") {
        return std::make_shared<PartialReplacement>(1.0, m_cache);  // Full replacement is just partial with ratio 1.0
    } else if (replacementMethod == "partial") {
        return std::make_shared<PartialReplacement>(replacementRatio, m_cache);
    }
    // Default partial
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
              const DNAInstance& instance,
              int generation,
              double bestFitness) const override {
        auto now = std::chrono::steady_clock::now();
        auto elapsedSec = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
        return (elapsedSec >= limitSeconds);
    }

    void reset() override {
        start = std::chrono::steady_clock::now();
    }

private:
    int limitSeconds;
    std::chrono::steady_clock::time_point start;
};

std::shared_ptr<IStopping> GAConfig::getStopping() const
{
    std::shared_lock<std::shared_mutex> lock(configMutex);
    
    // Get maxGenerations value safely under lock
    int maxGen = m_maxGenerations;
    std::cout << "[GAConfig] Creating stopping criteria with maxGenerations=" << maxGen << std::endl;
    
    if (stoppingMethod == "maxGenerations") {
        // Create with direct value instead of reference
        auto stopping = std::make_shared<MaxGenerationsStopping>(maxGen);
        std::cout << "[GAConfig] Created MaxGenerationsStopping using value: " << maxGen << std::endl;
        return stopping;
    } else if (stoppingMethod == "noImprovement") {
        return std::make_shared<NoImprovementStopping>(noImprovementGenerations);
    } else if (stoppingMethod == "timeLimit") {
        return std::make_shared<TimeLimitStopping>(timeLimitSeconds);
    }
    
    // Default case
    std::cout << "[GAConfig] Using default stopping criteria" << std::endl;
    auto stopping = std::make_shared<MaxGenerationsStopping>(maxGen);
    std::cout << "[GAConfig] Created default MaxGenerationsStopping using value: " << maxGen << std::endl;
    return stopping;
}

void GAConfig::setParameters(const ParameterSet& ps) {
    std::unique_lock<std::shared_mutex> lock(configMutex);
    
    // Set all parameters from the parameter set
    if (ps.contains("populationSize")) setPopulationSize(ps.getInt("populationSize"));
    if (ps.contains("mutationRate")) setMutationRate(ps.getDouble("mutationRate"));
    if (ps.contains("maxGenerations")) setMaxGenerations(ps.getInt("maxGenerations"));
    if (ps.contains("replacementRatio")) setReplacementRatio(ps.getDouble("replacementRatio"));
    if (ps.contains("crossoverProbability")) setCrossoverProbability(ps.getDouble("crossoverProbability"));
    if (ps.contains("tournamentSize")) setTournamentSize(ps.getInt("tournamentSize"));
    
    // Update string parameters
    if (ps.contains("selectionMethod")) selectionMethod = ps.getString("selectionMethod");
    if (ps.contains("crossoverType")) crossoverType = ps.getString("crossoverType");
    if (ps.contains("mutationMethod")) mutationMethod = ps.getString("mutationMethod");
    if (ps.contains("replacementMethod")) replacementMethod = ps.getString("replacementMethod");
    if (ps.contains("stoppingMethod")) stoppingMethod = ps.getString("stoppingMethod");
    if (ps.contains("fitnessType")) fitnessType = ps.getString("fitnessType");
}

bool GAConfig::validate() const {
    std::shared_lock<std::shared_mutex> lock(configMutex);
    
    // Check numeric ranges
    if (m_populationSize <= 0) return false;
    if (m_mutationRate < 0.0 || m_mutationRate > 1.0) return false;
    if (m_maxGenerations <= 0) return false;
    if (replacementRatio < 0.0 || replacementRatio > 1.0) return false;
    if (crossoverProbability < 0.0 || crossoverProbability > 1.0) return false;
    if (tournamentSize <= 0) return false;
    
    // Check string parameters
    if (selectionMethod.empty()) return false;
    if (crossoverType.empty()) return false;
    if (mutationMethod.empty()) return false;
    if (replacementMethod.empty()) return false;
    if (stoppingMethod.empty()) return false;
    if (fitnessType.empty()) return false;
    
    return true;
}

int GAConfig::getParentCount() const {
    return 2;  // Default to 2 parents for crossover
}

double GAConfig::getTargetFitness() const {
    return targetFitness;
}

int GAConfig::getTournamentSize() const {
    return tournamentSize;
}

