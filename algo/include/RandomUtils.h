#ifndef RANDOM_UTILS_H
#define RANDOM_UTILS_H

#include <random>

inline int random_node(int num_nodes) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, num_nodes - 1);
    return dis(gen);
}

#endif // RANDOM_UTILS_H
