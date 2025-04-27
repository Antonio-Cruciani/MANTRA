#ifndef TEMPORAL_DIAMETER_H
#define TEMPORAL_DIAMETER_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <cstdint> // for int64_t
#include <unordered_map>
#include <queue>
#include "TemporalGraph.h"
#include <tuple>

class TemporalDiameter {
public:
      // pointer to the temporal graph

    // Constructor
    TemporalDiameter(const TemporalGraph* tg, int64_t nn, int64_t ntn);

    // Declare methods you will implement (example placeholders)
    void run_traversal(double & diameter);

private:
    const TemporalGraph* tg;
    std::queue<std::tuple<int64_t, int64_t>> forward_queue;

    std::vector<int64_t> dist;     // distance per node
    std::vector<int64_t> dist_t;   // distance per temporal node

    // Any other private members you want later
};


#endif // TEMPORAL_DIAMETER_H
