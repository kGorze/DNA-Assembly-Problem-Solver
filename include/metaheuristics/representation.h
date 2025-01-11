#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include "generator/dna_generator.h"

enum class RepresentationType {
 DIRECT,      
 PERMUTATION, 
 GRAPH_PATH   
};

class IRepresentation {
public:
 virtual ~IRepresentation() = default;

 virtual std::vector<std::shared_ptr<std::vector<int>>> 
 initializePopulation(int popSize, const DNAInstance &instance) = 0;

 virtual std::string decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                                 const DNAInstance &instance) = 0;
};

class DirectDNARepresentation : public IRepresentation {
public:
 DirectDNARepresentation(int fallbackN = 500, double randomRatio = 1.0)
     : fallbackN(fallbackN), randomRatio(randomRatio) {}
  
 std::vector<std::shared_ptr<std::vector<int>>> 
 initializePopulation(int popSize, const DNAInstance &instance) override;

 std::string decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                         const DNAInstance &instance) override;

private:
 int fallbackN;
 double randomRatio;
};

class PermutationRepresentation : public IRepresentation {
public:
 std::vector<std::shared_ptr<std::vector<int>>> 
 initializePopulation(int popSize, const DNAInstance &instance) override;
    
 std::string decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                         const DNAInstance &instance) override;
};

class GraphPathRepresentation : public IRepresentation {
public:
 std::vector<std::shared_ptr<std::vector<int>>> 
 initializePopulation(int popSize, const DNAInstance &instance) override;
    
 std::string decodeToDNA(std::shared_ptr<std::vector<int>> individual, 
                         const DNAInstance &instance) override;
};

#endif //REPRESENTATION_H
