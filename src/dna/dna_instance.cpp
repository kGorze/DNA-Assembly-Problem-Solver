#include "dna/dna_instance.h"
#include "utils/logging.h"
#include <algorithm>

int DNAInstance::findStartVertexIndex(const DNAInstance& instance) {
    if (instance.getSpectrum().empty() || instance.getDNA().empty() || instance.getK() <= 0) {
        return -1;
    }

    std::string startFrag = instance.getDNA().substr(0, instance.getK());
    const auto& spectrum = instance.getSpectrum();
    
    for (int i = 0; i < static_cast<int>(spectrum.size()); i++) {
        if (spectrum[i] == startFrag) {
            return i;
        }
    }
    
    return -1;
} 