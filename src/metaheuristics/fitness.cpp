//
// Created by konrad_guest on 07/01/2025.
// SMART

#include <unordered_set>
#include <iostream>
#include "metaheuristics/fitness.h"

/* ================== SimpleFitness ================== */
double SimpleFitness::evaluate(std::shared_ptr<std::vector<int>> individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // decode
    std::string dna = representation->decodeToDNA(individual, instance);
    int k = instance.getK();
    if(k <= 0) return 0.0;

    // build set from dna's k-mers
    std::unordered_set<std::string> kmerSet;
    for(int i=0; i+(k) <= (int)dna.size(); i++){
        kmerSet.insert(dna.substr(i, k));
    }

    const auto &spec = instance.getSpectrum();
    std::unordered_set<std::string> spectrumSet(spec.begin(), spec.end());

    int matches = 0;
    for(const auto &ss : kmerSet) {
        if(spectrumSet.count(ss)) {
            matches++;
        }
    }

    return (double)matches;
}

/* ================== BetterFitness ================== */
double BetterFitness::evaluate(std::shared_ptr<std::vector<int>> individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    std::string dna = representation->decodeToDNA(individual, instance);

    int k = instance.getK();
    if (k <= 0) return 0.0;

    std::unordered_set<std::string> kmerSet;
    for(int i=0; i + k <= (int)dna.size(); i++){
        kmerSet.insert(dna.substr(i, k));
    }

    const auto &spec = instance.getSpectrum();
    std::unordered_set<std::string> spectrumSet(spec.begin(), spec.end());

    int matches = 0;
    for(auto &ss : kmerSet) {
        if(spectrumSet.count(ss)) {
            matches++;
        }
    }

    // Dodajmy wpływ odległości Levenshteina
    std::string originalDNA = instance.getDNA();
    if (originalDNA.empty()) {
        return (double)matches;
    }

    int distance = levenshteinDistance(originalDNA, dna);

    double alpha = 0.5; 
    double fitnessVal = matches - alpha * distance;
    return fitnessVal;
}

/* ================== SmithWatermanFitness ================== */
double SmithWatermanFitness::evaluate(std::shared_ptr<std::vector<int>> individual,
                                      const DNAInstance& instance,
                                      std::shared_ptr<IRepresentation> representation) const
{
    std::string reconstructed = representation->decodeToDNA(individual, instance);
    std::string original = instance.getDNA();
        
    if (reconstructed.empty() || original.empty()) {
        return 0.0;
    }
        
    int swScore = smithWaterman(original, reconstructed);
    double normalizedScore = swScore / 5.0; // np. dzielimy przez MATCH_SCORE
    return normalizedScore;
}

int SmithWatermanFitness::smithWaterman(const std::string& seq1, const std::string& seq2) const {
    int m = (int)seq1.length();
    int n = (int)seq2.length();
    
    std::vector<int> score(m+1, 0);
    int maxScore = 0;
    
    for (int j = 1; j <= n; j++) {
        int prev = 0;
        for (int i = 1; i <= m; i++) {
            int temp = score[i];
            int match = prev + (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
            int delete_gap = score[i] + GAP_PENALTY;
            int insert_gap = score[i-1] + GAP_PENALTY;
            
            score[i] = std::max({0, match, delete_gap, insert_gap});
            maxScore = std::max(maxScore, score[i]);
            prev = temp;
        }
    }
    return maxScore;
}

/* ================== OptimizedGraphBasedFitness ================== */

// Pomocnicza inicjalizacja bufora
void OptimizedGraphBasedFitness::initBuffers(size_t size) const {
    if (nodeUsageBuffer.size() < size) {
        nodeUsageBuffer.resize(size);
    }
}

// Generowanie klucza do kesza
std::string OptimizedGraphBasedFitness::createCacheKey(const std::vector<std::string>& spectrum, int k) const {
    std::string key = std::to_string(k) + ":";
    for (const auto& s : spectrum) {
        key += s + ",";
    }
    return key;
}

// Właściwe budowanie grafu (lub użycie z kesza)
std::vector<std::vector<OptimizedGraphBasedFitness::Edge>> 
OptimizedGraphBasedFitness::buildSpectrumGraph(
    const std::vector<std::string>& spectrum, 
    int k) const 
{
    std::string cacheKey = createCacheKey(spectrum, k);
    auto it = graphCache.find(cacheKey);
    if (it != graphCache.end()) {
        return it->second;
    }
    
    int n = (int)spectrum.size();
    std::vector<std::vector<Edge>> graph(n);
    
    // Równoległe wypełnianie list sąsiedztwa
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        graph[i].reserve(n/4); // np. heurystycznie
        for (int j = 0; j < n; j++) {
            if (i != j) {
                int weight = calculateEdgeWeight(spectrum[i], spectrum[j], k);
                if (weight > 0) {
                    // W sekcji równoległej dostęp do graph[i] jest bezpieczny,
                    // bo i jest unikalne w każdej iteracji i.
                    graph[i].emplace_back(j, weight);
                }
            }
        }
    }
    
    // Zapamiętujemy w keszu
    graphCache[cacheKey] = graph;
    return graph;
}

// Wyznaczenie wagi krawędzi między dwoma k-merami
int OptimizedGraphBasedFitness::calculateEdgeWeight(
    const std::string& from, 
    const std::string& to, 
    int k) const 
{
    // Podobna logika co wcześniej
    const char* fromEnd = from.data() + from.size() - (k - 1);
    const char* toStart = to.data();
    
    bool match = true;
    for (int i = 0; i < k - 1; i++) {
        if (fromEnd[i] != toStart[i]) {
            match = false;
            break;
        }
    }
    if (match) return 1;
    
    // sprawdzamy jeszcze "kolejne" przesunięcia ...
    fromEnd++;
    match = true;
    for (int i = 0; i < k - 2; i++) {
        if (fromEnd[i] != toStart[i]) {
            match = false;
            break;
        }
    }
    if (match) return 2;
    
    fromEnd++;
    // Jeszcze minimalna zbieżność – 1 znak
    if (fromEnd[0] == toStart[0]) return 3;
    
    return 0;
}

// Analiza ścieżki w grafie: w tej nowej wersji adjacencyMatrix budujemy lokalnie.
OptimizedGraphBasedFitness::PathAnalysis 
OptimizedGraphBasedFitness::analyzePath(const std::vector<int>& path,
                                        const std::vector<std::vector<PreprocessedEdge>>& adjacencyMatrix) const 
{
    // Wyzerowanie bufora
    std::fill(nodeUsageBuffer.begin(), nodeUsageBuffer.end(), 0);

    PathAnalysis analysis(nodeUsageBuffer);

    for (size_t i = 0; i < path.size() - 1; i++) {
        int from = path[i];
        int to = path[i + 1];
        
        if (from >= 0 && from < (int)adjacencyMatrix.size() &&
            to   >= 0 && to   < (int)adjacencyMatrix.size())
        {
            // zliczamy użycie węzła
            analysis.nodeUsageCount[from]++;
            if (analysis.nodeUsageCount[from] == 1) {
                analysis.uniqueNodesUsed++;
            } else {
                analysis.repeatNodeUsages++;
            }

            // sprawdzamy krawędź
            const auto& edge = adjacencyMatrix[from][to];
            if (edge.exists) {
                if (edge.weight == 1) {
                    analysis.edgesWeight1++;
                } else {
                    analysis.edgesWeight2or3++;
                }
            }
        }
    }
    // Ostatni węzeł
    if (!path.empty()) {
        int last = path.back();
        if (last >= 0 && last < (int)adjacencyMatrix.size()) {
            analysis.nodeUsageCount[last]++;
            if (analysis.nodeUsageCount[last] == 1) {
                analysis.uniqueNodesUsed++;
            } else {
                analysis.repeatNodeUsages++;
            }
        }
    }

    return analysis;
}

// Prosta funkcja – w permutacji 1:1 to nasza ścieżka.
const std::vector<int>& 
OptimizedGraphBasedFitness::permutationToPath(std::shared_ptr<std::vector<int>> individual) const {
    return *individual;
}

double OptimizedGraphBasedFitness::evaluate(std::shared_ptr<std::vector<int>> individual,
                                            const DNAInstance& instance,
                                            std::shared_ptr<IRepresentation> representation) const 
{
    const auto& spectrum = instance.getSpectrum();
    int k = instance.getK();
    if (spectrum.empty() || k <= 0) return 0.0;
    
    initBuffers(spectrum.size());
    // Budujemy (albo wczytujemy z kesza) listę sąsiedztwa:
    auto graph = buildSpectrumGraph(spectrum, k);

    // Aby analizować ścieżkę, potrzebujemy "macierzy" adjacencyMatrix, 
    // ale zrobimy ją tymczasowo w kontekście TYLKO TEGO wywołania evaluate().
    // Nie zapisujemy do pola globalnego klasy!
    int n = (int)graph.size();
    std::vector<std::vector<PreprocessedEdge>> adjacencyMatrix(n, std::vector<PreprocessedEdge>(n));

    // Lokalnie wypełniamy adjacencyMatrix
    // (to jest zamiast poprzedniego preprocessGraph(...) które modyfikowało pole klasy)
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; i++) {
        for (auto &edge : graph[i]) {
            adjacencyMatrix[i][edge.to] = PreprocessedEdge(edge.to, edge.weight, true);
        }
    }

    // Konwersja osobnika na ścieżkę
    const auto& path = permutationToPath(individual);

    // Analiza ścieżki
    auto analysis = analyzePath(path, adjacencyMatrix);
    
    // Wyliczenie finalnej oceny fitness
    double edgeScore = analysis.edgesWeight1 - (2.0 * analysis.edgesWeight2or3);
    double coverageScore = analysis.uniqueNodesUsed - (0.5 * analysis.repeatNodeUsages);
    
    const double alpha = 0.7;
    const double beta = 0.3;
    
    return (alpha * edgeScore) + (beta * coverageScore);
}
