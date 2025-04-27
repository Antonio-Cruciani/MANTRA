#include "TemporalGraph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;



TemporalGraph::TemporalGraph( const string &filename, const string& sep ) 
{
    load_from_file(filename,sep);
}


bool TemporalGraph::load_from_file(const string& filename, const string& sep) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "The temporal edge list file " << filename << " does not exist" << endl;
        return false;
    }

    int current_node_id = 0;
    int current_time = 0;
    unordered_map<string, int> file_id_to_graph_id;
    unordered_map<int, int> file_time_to_graph_time;
    set<tuple<int, int, int>> edge_set;

    string line;
    while (getline(infile, line)) {
        stringstream ss(line);
        string u, v, t_str;
        if (!getline(ss, u, sep[0]) || !getline(ss, v, sep[0]) || !getline(ss, t_str)) {
            cerr << "Bad line format: " << line << endl;
            continue;
        }

        int t = stoi(t_str);
        if (u != v) {
            if (file_id_to_graph_id.count(u) == 0) {
                file_id_to_graph_id[u] = current_node_id++;
                file_id.push_back(u);
            }
            if (file_id_to_graph_id.count(v) == 0) {
                file_id_to_graph_id[v] = current_node_id++;
                file_id.push_back(v);
            }
            if (file_time_to_graph_time.count(t) == 0) {
                file_time_to_graph_time[t] = current_time++;
                file_time.push_back(t);
            }
        }
    }

    infile.close();
    sort(file_time.begin(), file_time.end());
    for (size_t i = 0; i < file_time.size(); ++i) {
        file_time_to_graph_time[file_time[i]] = static_cast<int>(i);
    }

    infile.open(filename);
    while (getline(infile, line)) {
        stringstream ss(line);
        string u, v, t_str;
        getline(ss, u, sep[0]);
        getline(ss, v, sep[0]);
        getline(ss, t_str);

        if (u != v) {
            int t = stoi(t_str);
            edge_set.emplace(
                file_id_to_graph_id[u],
                file_id_to_graph_id[v],
                file_time_to_graph_time[t]
            );
        }
    }

    infile.close();
    temporal_edges.assign(edge_set.begin(), edge_set.end());
    sort(temporal_edges.begin(), temporal_edges.end(),
        [](const auto& a, const auto& b) { return get<2>(a) < get<2>(b); });

    num_nodes = static_cast<int>(file_id_to_graph_id.size());
    return true;
}

void TemporalGraph::print_stats(const string& graph_name) const {
    cout << "====================================================\n";
    cout << "Temporal network: " << graph_name << "\n";
    cout << "====================================================\n";
    cout << "Number of nodes: " << num_nodes << "\n";
    cout << "Number of temporal edges: " << temporal_edges.size() << "\n";
    cout << "Number of unique time stamps: " << file_time.size() << "\n";
    cout << "====================================================\n";
}

vector<vector<pair<int, int>>> TemporalGraph::temporal_adjacency_list() const {
    vector<vector<pair<int, int>>> adjacency_list(num_nodes);
    for (const auto& edge : temporal_edges) {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int t = get<2>(edge);
        adjacency_list[u].emplace_back(v, t);
    }
    return adjacency_list;
}

vector<vector<pair<int, int>>> TemporalGraph::temporal_incidency_list() const {
    vector<vector<pair<int, int>>> incidence_list(num_nodes);
    for (const auto& edge : temporal_edges) {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int t = get<2>(edge);
        incidence_list[v].emplace_back(u, t);
    }
    return incidence_list;
}

vector<pair<int, int>> TemporalGraph::next_temporal_neighbors(
    const vector<vector<pair<int, int>>>& tal,
    int u,
    int t
) const {
    const vector<pair<int, int>>& neighbors = tal[u];
    int left = 0;
    int right = neighbors.size() - 1;
    int pos = neighbors.size();

    while (left <= right) {
        int mid = (left + right) / 2;
        if (neighbors[mid].second <= t) {
            left = mid + 1;
        } else {
            pos = mid;
            right = mid - 1;
        }
    }

    return vector<pair<int, int>>(neighbors.begin() + pos, neighbors.end());
}

unordered_map<pair<int, int>, int, hash<pair<int, int>>> TemporalGraph::temporal_node_index() const {
    unordered_map<pair<int, int>, int, hash<pair<int, int>>> index_map;
    int current_index = 0;

    for (const auto& edge : temporal_edges) {
        //int u = get<0>(edge);
        int v = get<1>(edge);
        int t = get<2>(edge);
        //pair<int, int> key_u = make_pair(u, t);
        pair<int, int> key_v = make_pair(v, t);

        //if (index_map.find(key_u) == index_map.end()) {
        //    index_map[key_u] = current_index++;
       // }
        if (index_map.find(key_v) == index_map.end()) {
            index_map[key_v] = current_index++;
        }
    }

    for (int s = 0; s < num_nodes; ++s) {
        pair<int, int> zero_key = make_pair(s, 0);
        if (index_map.find(zero_key) == index_map.end()) {
            index_map[zero_key] = current_index++;
        }
    }
    return index_map;
}

