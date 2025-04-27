#include "TemporalDiameter.h"
#include <iostream>  // only if you plan to print debug info
#include "Sampler.h"
#include <queue>
#include "RandomUtils.h"
using namespace std;
TemporalDiameter::TemporalDiameter(const TemporalGraph* tg, int64_t nn, int64_t ntn)
    : tg(tg),            // store the pointer to the graph
      dist(nn, -1),        // initialize dist with nn entries, all -1 (meaning unvisited)
      dist_t(ntn, -1)      // initialize dist_t with ntn entries, all -1
{
    
    // Constructor body (optional debug info)
    // std::cout << "TemporalDiameter initialized with " << nn << " nodes and " << ntn << " temporal nodes.\n";
}




void TemporalDiameter::run_traversal(double& diameter) {
    auto tal = tg->temporal_adjacency_list();
    auto tn_index = tg->temporal_node_index();
    int s,z;
    unordered_set<int> active_times_for_z;
    s = random_node(tg->num_nodes);
  
    fill(dist.begin(), dist.end(), -1);
    fill(dist_t.begin(), dist_t.end(), -1);
    while (!forward_queue.empty()) forward_queue.pop();
    // Initialize BFS
    int64_t tni = tn_index[make_pair(s, 0)];
    dist[s] = 0;
    dist_t[tni] = 0;
    forward_queue.emplace(s, 0);
    double max_distance = 0;
    
    
    while (!forward_queue.empty()) {
        auto [u, t] = forward_queue.front(); forward_queue.pop();
        int64_t tni_u = tn_index[make_pair(u, t)];


        for (const auto& [w, t_w] : tg->next_temporal_neighbors(tal, u, t)) {
            int64_t tni_w = tn_index[make_pair(w, t_w)];

            if (dist_t[tni_w] == -1) {
                dist_t[tni_w] = dist_t[tni_u] + 1;
                if (dist[w] == -1) {
                    dist[w] = dist_t[tni_u] + 1;
                    if ( dist[w]> max_distance){
                        max_distance =  dist[w];
                    }
                }
                forward_queue.emplace(w, t_w);
            }

           
        }

    }

    #pragma omp critical
    {
        if (max_distance> diameter){
            diameter = max_distance;
        }
       
    }
    return;
}
    