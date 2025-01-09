//
// Created by konrad_guest on 09/01/2025.
//

#ifndef ZOBRIST_HASHER_H
#define ZOBRIST_HASHER_H

#include <vector>
#include <random>
#include <cstdint>

class ZobristHasher {
public:
    static ZobristHasher& getInstance() {
        static ZobristHasher instance;
        return instance;
    }
    
    // Hash dla permutacji (główny typ reprezentacji w obecnym kodzie)
    uint64_t hashPermutation(const std::vector<int>& perm) {
        uint64_t hash = 0;
        for (size_t i = 0; i < perm.size(); ++i) {
            hash ^= randomTable[i][perm[i]];
        }
        return hash;
    }
    
private:
    ZobristHasher() {
        initializeRandomTable();
    }
    
    void initializeRandomTable() {
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint64_t> dis;
        
        // Inicjalizacja tablicy dla możliwych pozycji i wartości
        // Zakładamy maksymalną długość permutacji 1000 i wartości 0-999
        randomTable.resize(1000);
        for (auto& row : randomTable) {
            row.resize(1000);
            for (auto& val : row) {
                val = dis(gen);
            }
        }
    }
    
    std::vector<std::vector<uint64_t>> randomTable;
};

#endif //ZOBRIST_HASHER_H
