#ifndef PREPROCESSED_EDGE_H
#define PREPROCESSED_EDGE_H

struct PreprocessedEdge {
    int weight;
    int overlap;
    bool valid;
    
    PreprocessedEdge(int w = 0, int o = 0, bool v = false) 
        : weight(w), overlap(o), valid(v) {}
};

#endif // PREPROCESSED_EDGE_H 