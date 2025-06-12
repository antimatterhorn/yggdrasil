// Copyright (C) 2025  Cody Raskin

#include "eosTable.hh"

EOSTable::EOSTable(std::vector<std::vector<double>> table,
            std::vector<double> x, std::vector<double> y)
        : values(std::move(table)), xgrid(std::move(x)), ygrid(std::move(y)) {}