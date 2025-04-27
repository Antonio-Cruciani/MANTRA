#include <random>
#include "TemporalGraph.h"
#include "Sampler.h"
#include <unordered_map>
#include <unordered_set>

#include <iostream>

Sampler::Sampler(const TemporalGraph *tg ,int64_t nn, int64_t ntn): 
    tg(tg),
    sigma(nn),
    dist(nn),
    sigma_t(ntn),
    sigma_z(ntn, 0),
    dist_t(ntn),
    predecessors(ntn),
    boolean_matrix(ntn, false)
{


    // queues are default-initialized and empty
}
// Returns a random node.
// Returns a random node ID in [0, num_nodes - 1]




void Sampler::run_sample(vector<double>& tilde_b,vector<double>& sp_lengths) {
    auto tal = tg->temporal_adjacency_list();
    auto tn_index = tg->temporal_node_index();
    int s,z;
    unordered_set<int> active_times_for_z;
    int n = tg->num_nodes;
    s = random_node(n);
    z = random_node(n);
    while (s==z){
        z = random_node(n);
    }
    fill(sigma.begin(), sigma.end(), 0);
    fill(dist.begin(), dist.end(), -1);
    fill(sigma_t.begin(), sigma_t.end(), 0);
    fill(dist_t.begin(), dist_t.end(), -1);
    fill(boolean_matrix.begin(), boolean_matrix.end(), false);
    for (auto& pred : predecessors) pred.clear();
    while (!forward_queue.empty()) forward_queue.pop();
    while (!backward_queue.empty()) backward_queue.pop();
    // Initialize BFS
    int64_t tni = tn_index[make_pair(s, 0)];
    sigma[s] = 1;
    sigma_t[tni] = 1;
    dist[s] = 0;
    dist_t[tni] = 0;
    forward_queue.emplace(s, 0);
    int64_t d_z_min = INT64_MAX;
    
    
    while (!forward_queue.empty()) {
        auto [u, t] = forward_queue.front(); forward_queue.pop();
        int64_t tni_u = tn_index[make_pair(u, t)];

        if (dist_t[tni_u] >= d_z_min) continue;

        for (const auto& [w, t_w] : tg->next_temporal_neighbors(tal, u, t)) {
            int64_t tni_w = tn_index[make_pair(w, t_w)];

            if (dist_t[tni_w] == -1) {
                dist_t[tni_w] = dist_t[tni_u] + 1;
                if (dist[w] == -1) {
                    dist[w] = dist_t[tni_u] + 1;
                    if (w == z){ 
                        d_z_min = dist[w];
                        active_times_for_z.insert(t_w);
                    }
                }
                forward_queue.emplace(w, t_w);
            }

            if (dist_t[tni_w] == dist_t[tni_u] + 1) {
                sigma_t[tni_w] += sigma_t[tni_u];
                predecessors[tni_w].emplace(u, t);
                if (dist_t[tni_w] == dist[w]) {
                    sigma[w] += sigma_t[tni_u];
                }
            }
        }
    }
    
    //Backtracking
    if (dist[z] > 0) {
        #pragma omp critical
        {
            sp_lengths[dist[z]] += 1;
        }
        //std::cout<<"DIST > 0 "<<s<<" , "<<z<<std::endl;
        // === Backward Phase ===
        

        for (auto& val : sigma_z) val = 0;
        fill(boolean_matrix.begin(), boolean_matrix.end(), false);

        tni = tn_index[make_pair(s, 0)];
        sigma_z[tni] = 1;
        //std::cout<<"CIAO "<<std::endl;
        //size_t T = tg.file_time.size();
        //for (int64_t t = 0; t < static_cast<int64_t>(T); ++t) {
        for (int t : active_times_for_z) {
            auto it = tn_index.find(make_pair(z, t));
            if (it == tn_index.end()) continue;
            
            int idx = it->second;
            if (sigma_t[idx] > 0) {
                for (auto& pred : predecessors[idx]) {
                    int64_t tni_w = tn_index.at(make_pair(get<0>(pred), get<1>(pred)));
                    sigma_z[tni_w] += 1;
                    if (!boolean_matrix[tni_w]) {
                        backward_queue.push(pred);
                        boolean_matrix[tni_w] = true;
                    }
                }
            }
        }



        while (!backward_queue.empty()) {
            auto [u, t] = backward_queue.front(); backward_queue.pop();
            tni = tn_index[make_pair(u, t)];
            if (u != s) {
                #pragma omp critical
                {
                    tilde_b[u] += static_cast<double>(sigma_z[tni]) *
                                  (static_cast<double>(sigma_t[tni]) / static_cast<double>(sigma[z]));
                    //std::cout<<"INCREASING "<<u<<" of "<<tilde_b[u]<<endl;
                }
            
                for (auto& pred : predecessors[tni]) {
                    int64_t tni_w = tn_index[make_pair(get<0>(pred), get<1>(pred))];
                    sigma_z[tni_w] += sigma_z[tni];
                    if (!boolean_matrix[tni_w]) {
                        backward_queue.push(pred);
                        boolean_matrix[tni_w] = true;
                    }
                }
            }
        }
        
    }

    return;
}








string Sampler::uint128_to_string(__uint128_t value) {
    if (value == 0) return "0";

    string result;
    while (value > 0) {
        result = static_cast<char>('0' + value % 10) + result;
        value /= 10;
    }
    return result;
}