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

size_t findIndex(const std::vector<double>& grid, double value) {
    auto it = std::upper_bound(grid.begin(), grid.end(), value);
    if (it == grid.begin() || it == grid.end()) return std::max<size_t>(grid.size() - 2, 0);
    return std::distance(grid.begin(), it) - 1;
}

double 
bilinearInterp(const std::vector<std::vector<double>>& table, 
                const std::vector<double>& xgrid, const std::vector<double>& ygrid,
                double x, double y) {
    size_t i = Lin::findIndex(xgrid, x);
    size_t j = Lin::findIndex(ygrid, y);

    double x0 = xgrid[i], x1 = xgrid[i + 1];
    double y0 = ygrid[j],   y1 = ygrid[j + 1];

    double q11 = table[i][j];
    double q12 = table[i][j + 1];
    double q21 = table[i + 1][j];
    double q22 = table[i + 1][j + 1];

    double tx = (x - x0) / (x1 - x0);
    double ty = (y - y0) / (y1 - y0);

    return (1 - tx) * (1 - ty) * q11 +
            (1 - tx) * ty       * q12 +
            tx       * (1 - ty) * q21 +
            tx       * ty       * q22;
}


}