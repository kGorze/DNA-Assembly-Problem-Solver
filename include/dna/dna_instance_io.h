#pragma once

#include "dna_instance.h"
#include <string>
#include <mutex>
#include <stdexcept>
#include <fstream>
#include <sstream>

class InstanceIO {
public:
    static bool saveInstance(const DNAInstance& instance, const std::string& filename);
    static DNAInstance loadInstance(const std::string& filename);
}; 