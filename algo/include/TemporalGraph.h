#ifndef TEMPORAL_GRAPH_H
#define TEMPORAL_GRAPH_H

#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <unordered_map>
#include <utility>

class TemporalGraph {
public:

    int num_nodes;
    std::vector<std::tuple<int, int, int>> temporal_edges;
    std::vector<std::string> file_id;
    std::vector<int> file_time;
    TemporalGraph( const std::string& filename, const std::string& sep);
    bool load_from_file(const std::string& filename, const std::string& sep = " ");
    void print_stats(const std::string& graph_name = "anonymous") const;

    std::vector<std::vector<std::pair<int, int>>> temporal_adjacency_list() const;
    std::vector<std::vector<std::pair<int, int>>> temporal_incidency_list() const;
    std::vector<std::pair<int, int>> next_temporal_neighbors(
        const std::vector<std::vector<std::pair<int, int>>>& tal,
        int u, int t
    ) const;
    std::unordered_map<int, int> file_time_to_graph_time;

    std::unordered_map<std::pair<int, int>, int, std::hash<std::pair<int, int>>> temporal_node_index() const;
    int get_nn(){
        return num_nodes;
    }
    int get_tni(){
        return num_nodes;
    }
};

// --- Hash specialization with guard ---
#ifndef HASH_PAIR_INT_INT_DEFINED
#define HASH_PAIR_INT_INT_DEFINED
namespace std {
    template <>
    struct hash<std::pair<int, int>> {
        size_t operator()(const std::pair<int, int>& p) const {
            return hash<int>()(p.first) ^ (hash<int>()(p.second) << 1);
        }
    };
}
#endif

#endif // TEMPORAL_GRAPH_H
