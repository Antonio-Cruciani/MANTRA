#ifndef SAMPLER_H
#define SAMPLER_H

#include <vector>
#include <set>
#include <tuple>
#include <queue>
#include <cstdint>
#include "TemporalGraph.h"  
#include "RandomUtils.h"
using namespace std;

class Sampler {
public:
    vector<__uint128_t> sigma;
    vector<int64_t> dist;
    vector<__uint128_t> sigma_t;
    vector<__uint128_t> sigma_z;
    vector<int64_t> dist_t;
    vector<set<tuple<int64_t, int64_t>>> predecessors;

    vector<bool> boolean_matrix;
    queue<tuple<int64_t, int64_t>> forward_queue;
    queue<tuple<int64_t, int64_t>> backward_queue;

    // Constructor
    Sampler(const TemporalGraph* tg, int64_t nn, int64_t ntn);

    // Utility
    //inline int random_node(int num_nodes);

    // Core logic
    void run_sample( vector<double>& tilde_b, vector<double>& sp_lengths);
    string uint128_to_string(__uint128_t value);
private:
    const TemporalGraph* tg;


};

#endif // SAMPLER_H
