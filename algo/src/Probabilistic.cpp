#include <random>
#include "utilities.h"  
#include "TemporalGraph.h"
#include "Sampler.h"
#include "TemporalDiameter.h"
#include <unordered_map>
#include <unordered_set>
#include "Probabilistic.h"
#include <iostream>
#include <omp.h>
#include <time.h>
#include <fstream>
#include <algorithm>
#include <cmath> 
Probabilistic::Probabilistic(const std::string& filename, const std::string& sep)
: tg(filename, sep)
{
    tilde_b.resize(tg.get_nn(), 0.0);
}

void Probabilistic::run_algorithm(){
    int n = tg.get_nn();
    //vector<double> tilde_b(n, 0.0);
    auto index_map = tg.temporal_node_index();
    int sample_size = 1000;
    // Create a single Sampler object
    int64_t max_index = 0;
    for (const auto& [_, idx] : index_map) {
        max_index = std::max(max_index, static_cast<int64_t>(idx));
    }
    cout<<"INDEX SIZE "<<index_map.size()<<" max index +1 "<<max_index+1<<endl;
    Sampler sampler(&tg,n, max_index+1);
    cout<<"MIAOO"<<endl;
    cout<<"SAMPLE SIZE "<<sample_size<<endl;
    // === Run samples in parallel ===
    double diameter = 0;
    omp_set_num_threads(20);
    double sample_size_diam = 256;
    //double start_time_diam = get_time_sec();
    #pragma omp parallel
    {
      for (int i = 0; i < sample_size_diam; ++i) {
        TemporalDiameter diam(&tg, n,max_index+1);
        diam.run_traversal(diameter);
      }
    }
    //double finish_time_diam = get_time_sec() - start_time_diam;
    //cout<<"Estimated Temporal Diameter "<<diameter<<" in "<<finish_time_diam<<" seconds"<<endl;
    cout<<"Estimated Temporal Diameter "<<diameter<<endl;

    // Bootstrap phase
    vector<double> sp_lengths(n,0.0);
    //double start_time_bp = get_time_sec();
    #pragma omp parallel
      {
        for (int i = 0; i < sample_size; ++i) {
            Sampler sampler(&tg, n,max_index+1);
            sampler.run_sample(tilde_b,sp_lengths);
        }
    }
    //double finish_time_bp = get_time_sec() - start_time_bp;
    //cout<<"Bootstrap Phase Completed in "<<finish_time_bp<<" seconds"<<endl;
     double avg_diam_ = 0.;
    for( int i=0; i <= n; i++ ){
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
    for( int i=0; i <= n; i++ ){
      for( int j=i+1; j <= n; j++ ){
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
      for (int u = 0; u < n; ++u) {
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

      return;

}


void Probabilistic::save_tilde_b(const vector<double>& tilde_b, const string& filename) {
  std::ofstream out(filename);
  if (!out.is_open()) {
      std::cerr << "Failed to open file " << filename << " for writing." << std::endl;
      return;
  }

  for (size_t i = 0; i < tilde_b.size(); ++i) {
      out << tilde_b[i] << "\n";
  }

  out.close();
  std::cout << "Saved tilde_b to " << filename << std::endl;
}