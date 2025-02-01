#ifndef OPTYMALIZACJA_KOMBINATORYCZNA_INSTANCE_IO_H
#define OPTYMALIZACJA_KOMBINATORYCZNA_INSTANCE_IO_H

#include "dna/dna_instance.h"
#include <string>

class InstanceIO {
public:
    virtual ~InstanceIO() = default;
    
    // Save DNAInstance to a file
    virtual bool saveInstance(const DNAInstance& instance, const std::string& filename) = 0;
    
    // Load DNAInstance from a file
    virtual DNAInstance loadInstance(const std::string& filename) = 0;
};

#endif //OPTYMALIZACJA_KOMBINATORYCZNA_INSTANCE_IO_H 