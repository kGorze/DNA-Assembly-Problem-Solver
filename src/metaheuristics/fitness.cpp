//
// Created by konrad_guest on 07/01/2025.
//



// ================== SimpleFitness ==================
#include <unordered_set>
#include <iostream>

#include "metaheuristics/fitness.h"


double SimpleFitness::evaluate(void* individual,
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

    // compare with original instance's spectrum
    const auto &spec = instance.getSpectrum();
    std::unordered_set<std::string> spectrumSet(spec.begin(), spec.end());

    // count how many from kmerSet are in spectrum
    int matches = 0;
    for(auto &ss : kmerSet) {
        if(spectrumSet.count(ss)) {
            matches++;
        }
    }

    return (double)matches;
}

double BetterFitness::evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // decode
    std::string dna = representation->decodeToDNA(individual, instance);

    // count matched k-mers
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

    // let's incorporate Levenshtein
    std::string originalDNA = instance.getDNA();
    if (originalDNA.empty()) {
        return (double)matches;
    }

    int distance = levenshteinDistance(originalDNA, dna);

    // synergy measure
    double alpha = 0.5; // adjust as needed
    double fitnessVal = matches - alpha * distance;
    return fitnessVal;
}

double SmithWatermanFitness::evaluate(void* individual,
                               const DNAInstance &instance,
                               std::shared_ptr<IRepresentation> representation) const
{
    // Get the DNA sequences
    std::string reconstructed = representation->decodeToDNA(individual, instance);
    std::string original = instance.getDNA();
        
    if (reconstructed.empty() || original.empty()) {
        return 0.0;
    }
        
    // Calculate Smith-Waterman score
    int swScore = smithWaterman(original, reconstructed);
        
    // Divide by 5 as mentioned in the paper for comparison purposes
    double normalizedScore = swScore / 5.0;
        
    return normalizedScore;
}

int SmithWatermanFitness::smithWaterman(const std::string& seq1, const std::string& seq2) const {
    int m = seq1.length();
    int n = seq2.length();
        
    // Initialize scoring matrix
    std::vector<std::vector<int>> score(m + 1, std::vector<int>(n + 1, 0));
        
    // Fill the scoring matrix
    int maxScore = 0;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            // Calculate match/mismatch score
            int match = score[i-1][j-1] + 
                (seq1[i-1] == seq2[j-1] ? MATCH_SCORE : MISMATCH_SCORE);
                
            // Calculate gap scores
            int delete_gap = score[i-1][j] + GAP_PENALTY;
            int insert_gap = score[i][j-1] + GAP_PENALTY;
                
            // Take the maximum score
            score[i][j] = std::max({0, match, delete_gap, insert_gap});
                
            // Update maximum score seen so far
            maxScore = std::max(maxScore, score[i][j]);
        }
    }
        
    return maxScore;
}
