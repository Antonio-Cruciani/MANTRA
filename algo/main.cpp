#include <iostream>
#include <algorithm>
#include "Sampler.h"
#include "TemporalGraph.h"
#include <omp.h>
#include <fstream>
#include <cmath>
#include "TemporalDiameter.h"
#include "Probabilistic.h"
#include <string>
using namespace std;

// you need to use the other temporal index for s an z!!
int main() {
    string filename = "example.txt"; 
    string sep = " ";
    Probabilistic algo(filename,sep);
    algo.run_algorithm();
   
}

/*
void test_temporal_node_index() {
    TemporalGraph tg;
    tg.load_from_file("test_graph.txt");

    auto index_map = tg.temporal_node_index();

    std::cout << "Testing temporal_node_index...\n";

    // Expected keys
    std::vector<std::pair<int, int>> expected_keys;

    for (const auto& edge : tg.temporal_edges) {
        expected_keys.emplace_back(std::get<1>(edge), std::get<2>(edge)); // (v, t)
    }

    for (int s = 0; s < tg.num_nodes; ++s) {
        expected_keys.emplace_back(s, 0); // (s, 0)
    }

    for (const auto& key : expected_keys) {
        assert(index_map.find(key) != index_map.end() && "Key not found in index_map");
    }

    std::cout << "All expected keys are present in the index.\n";

    // Optional: Check uniqueness
    std::set<int> seen_indices;
    for (const auto& [key, idx] : index_map) {
        assert(seen_indices.insert(idx).second && "Duplicate index detected!");
    }

    std::cout << "All indices are unique.\n";
}
int main() {
    test_temporal_node_index();
    return 0;
}*/


/* string filename = "example.txt";  // <-- Replace with your file path
    if (!tg.load_from_file(filename)) {
        cerr << "Failed to load graph." << endl;
        return 1;
    }

    tg.print_stats("TestGraph");

    auto adj_list = tg.temporal_adjacency_list();

    // Print first 3 nodes' adjacency lists (if available)
    cout << "\nSample Temporal Adjacency List (first 3 nodes):\n";
    for (int u = 0; u < min(3, tg.num_nodes); ++u) {
        cout << "Node " << u << ":";
        for (const auto& [v, t] : adj_list[u]) {
            cout << " -> (" << v << ", t=" << t << ")";
        }
        cout << endl;
    }

    // Test next_temporal_neighbors for node 0 after time 2
    int test_node = 0;
    int test_time = 2;
    auto neighbors = tg.next_temporal_neighbors(adj_list, test_node, test_time);
    cout << "\nNeighbors of node " << test_node << " after time " << test_time << ":\n";
    for (const auto& [v, t] : neighbors) {
        cout << " -> (" << v << ", t=" << t << ")\n";
    }

    // Test temporal_node_index
    auto index_map = tg.temporal_node_index();
    cout << "\nTotal temporal node-time pairs indexed: " << index_map.size() << endl;

    // === Test next_temporal_neighbors for a few test cases
    vector<tuple<int, int>> test_queries = {
        {0, 0}, {0, 2}, {1, 1},{2, 1},{2, 2}, {2, 3},{2, 4},{2, 100}
    };

    cout << "\n=== Testing next_temporal_neighbors ===\n";
    for (const auto& [u, t] : test_queries) {
        cout << "Neighbors of node " << u << " after time " << t << ":\n";
        auto neighbors = tg.next_temporal_neighbors(adj_list, u, t);
        if (neighbors.empty()) {
            cout << "  (none)\n";
        } else {
            for (const auto& [v, time] : neighbors) {
                cout << "  -> (" << v << ", t=" << time << ")\n";
            }
        }
    }



    
    tg.print_stats("SampleGraph");

    int sample_size = 100;
    vector<double> tilde_b(tg.num_nodes, 0.0);  // Shared betweenness vector
    // Create a single Sampler object
    cout<<"bau"<<endl;
    int64_t max_index = 0;
    for (const auto& [_, idx] : index_map) {
        max_index = std::max(max_index, static_cast<int64_t>(idx));
    }
    cout<<"INDEX SIZE "<<tg.temporal_node_index().size()<<" max index +1 "<<max_index+1<<endl;
    Sampler sampler(&tg, tg.num_nodes, max_index+1);
    cout<<"MIAOO"<<endl;
    cout<<"SAMPLE SIZE "<<sample_size<<endl;
    // === Run samples in parallel ===
    vector<double> sp_lengths(tg.num_nodes,0.0);
    double diameter = 0;
    omp_set_num_threads(1);
    double sample_size_diam = 256;
    #pragma omp parallel
    {
      for (int i = 0; i < sample_size; ++i) {
        TemporalDiameter diam(&tg, tg.num_nodes,max_index+1);
        diam.run_traversal(tg, diameter);
      }
    }
    cout<<"Estimated Temporal Diameter "<<diameter<<endl;


    #pragma omp parallel
      {
        for (int i = 0; i < sample_size; ++i) {
            Sampler sampler(&tg, tg.num_nodes,max_index+1);
            sampler.run_sample(tg, tilde_b,sp_lengths);
        }
    }

     double avg_diam_ = 0.;
    for( int i=0; i <= tg.num_nodes; i++ ){
        avg_diam_ += (double)i*sp_lengths[i]/(double)(sample_size);
    }
    // upper bound using bernstein
    double delta = 0.05;
  double log_term_avgspl = log(1./delta);
  double const_term_avgspl = (diameter-1)*log_term_avgspl/(double)sample_size;
  double avg_diam_upperbound_b = avg_diam_ + const_term_avgspl + sqrt( 2*const_term_avgspl*avg_diam_ + pow(const_term_avgspl,2.) );
  bool debug_sp_lengths = false;

  // upper bound using empirical bernstein
  double var_estimate_avg_diam = 0.;
  for( int i=0; i <= tg.num_nodes; i++ ){
    for( int j=i+1; j <= tg.num_nodes; j++ ){
      var_estimate_avg_diam += pow(i-j,2.0)*sp_lengths[i]/(double)sample_size*sp_lengths[j]/(double)(sample_size-1);
    }
   
  }
  log_term_avgspl = log(2./delta);
  double avg_diam_upperbound_eb = avg_diam_ + 7./3.*(diameter-1)*log_term_avgspl/(double)sample_size + sqrt( 2*var_estimate_avg_diam*log_term_avgspl/(double)sample_size );

  double avg_diam_upperbound = min(avg_diam_upperbound_b , avg_diam_upperbound_eb);
  avg_diam_upperbound = min(avg_diam_upperbound , (double)(diameter-1));

    cout << "avg_diam_: " << avg_diam_ << std::endl;
    cout << "avg_diam_upperbound_b: " << avg_diam_upperbound_b << std::endl;
    cout << "avg_diam_upperbound_eb: " << avg_diam_upperbound_eb << std::endl;
    cout << "var_estimate_avg_diam: " << var_estimate_avg_diam << std::endl;
  
    // === Print top-k scores ===
    int k = 10;
    vector<pair<int, double>> centrality;
    for (int u = 0; u < tg.num_nodes; ++u) {
        centrality.emplace_back(u, tilde_b[u]);
    }

    // Sort by centrality descending
    sort(centrality.begin(), centrality.end(),
         [](const auto& a, const auto& b) { return a.second > b.second; });
    cout<<"LEN CENT SIZE "<<centrality.size()<<endl;
    cout << "\nTop " << k << " nodes by temporal centrality:\n";
    for (int i = 0; i < min(k, (int)centrality.size()); ++i) {
        cout << "Node " << centrality[i].first << ": " << centrality[i].second/sample_size << "\n";
    }
    for (int i = 0; i < tilde_b.size(); ++i) {
        tilde_b[i] = tilde_b[i]/sample_size;
    }
    cout<<"DONE"<<endl;
    save_tilde_b(tilde_b,"results_tilde_b.txt");
    return 0;
}
void save_tilde_b(const std::vector<double>& tilde_b, const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file " << filename << " for writing." << std::endl;
        return;
    }

    for (size_t i = 0; i < tilde_b.size(); ++i) {
        out << tilde_b[i] << "\n";
    }

    out.close();
    std::cout << "Saved tilde_b to " << filename << std::endl;*/