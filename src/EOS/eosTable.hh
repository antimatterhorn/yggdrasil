// Copyright (C) 2025  Cody Raskin

#pragma once

#include <vector>
#include "../Math/scalarMath.hh"

class EOSTable {
    std::vector<std::vector<double>> values;
    std::vector<double> xgrid, ygrid;

public:
    EOSTable() = default;
    EOSTable(std::vector<std::vector<double>> table,
            std::vector<double> x, std::vector<double> y);

    double interpolate(double x, double y) const {
        return Lin::bilinearInterp(values, xgrid, ygrid, x, y);
    }
};

#include "eosTable.cc"
