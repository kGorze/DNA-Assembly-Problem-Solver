//
// Created by konrad_guest on 23/01/2025.
//

#ifndef PARAMETERS_PARSER_H
#define PARAMETERS_PARSER_H

#include "racing.h"
#include "tuning_structures.h"
#include "meta_ea.h"
#include "parameter_tuning_manager.h"

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <iostream>

// Zakładam, że masz zdefiniowane struktury ParameterSet i TuningResult

class ParameterParser {
public:
    // Funkcje publiczne
    static ParameterSet parseConfigFile(const std::string &filename);
    
    static std::vector<ParameterSet> generateRandomCandidates(int numCandidates);
    static std::vector<ParameterSet> generateGridOfCandidatesWithout();
    static std::vector<ParameterSet> generateGridOfCandidatesRS();
};

#endif //PARAMETERS_PARSER_H
