#ifndef DNA_INSTANCE_IO_H
#define DNA_INSTANCE_IO_H

#include "dna_instance.h"
#include <string>

class InstanceIO {
public:
    static bool loadInstance(const std::string& filename, DNAInstance& instance);
    static bool saveInstance(const DNAInstance& instance, const std::string& filename);
};

#endif // DNA_INSTANCE_IO_H 