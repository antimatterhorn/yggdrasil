#pragma once

#include <vector>
#include <cmath>

namespace Lin {

std::vector<double> linspace(double Start, double End, size_t num) {
    std::vector<double> result(num);
    double delta = (End - Start) / (num - 1);
    for (size_t i = 0; i < num; ++i) {
        result[i] = Start + i * delta;
    }
    return result;
}

}