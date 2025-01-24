#ifndef PREPROCESSED_EDGE_H
#define PREPROCESSED_EDGE_H

struct PreprocessedEdge {
    int to;
    int weight;
    bool exists;
    PreprocessedEdge(int t = 0, int w = 0, bool e = false) 
        : to(t), weight(w), exists(e) {}
};

#endif // PREPROCESSED_EDGE_H 