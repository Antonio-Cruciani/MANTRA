#ifndef PROBABILISTIC_H
#define PROBABILISTIC_H

#include <string>
#include <vector>
#include <tuple>
#include <set>
#include <unordered_map>
#include <utility>
#include "utilities.h"
#include "TemporalGraph.h"

class Probabilistic{
    public:
        std::vector<double> tilde_b;
        // Constructor
        Probabilistic(const std::string & filename ,const std::string& sep);
        void run_algorithm();
        void save_tilde_b(const std::vector<double>& tilde_b, const std::string& filename);
    private:
        TemporalGraph tg;
};



#endif // PROBABILISTIC_H
