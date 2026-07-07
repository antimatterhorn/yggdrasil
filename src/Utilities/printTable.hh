// Copyright (C) 2026  Cody Raskin
#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "../DataBase/field.hh"

namespace PrintTableDetail {
    inline void collect(std::vector<const Field<double>*>&) {}

    template <typename... Rest>
    inline void collect(std::vector<const Field<double>*>& columns, const Field<double>* field, Rest... rest) {
        columns.push_back(field);
        collect(columns, rest...);
    }
}

// Prints a table of scalar Field<double> columns, one row per node/zone.
// Column headers come from each field's own name, e.g.:
//   printTable(rho->size(), r, rho, m, u, P, T);
template <typename... Fields>
void printTable(unsigned int numRows, const Field<double>* first, Fields... rest) {
    std::vector<const Field<double>*> columns;
    columns.push_back(first);
    PrintTableDetail::collect(columns, rest...);

    const int colWidth = 16;
    for (const auto* col : columns)
        std::cout << std::setw(colWidth) << (col->hasName() ? col->getNameString() : std::string("-"));
    std::cout << "\n";

    std::cout << std::scientific << std::setprecision(6);
    for (unsigned int i = 0; i < numRows; ++i) {
        for (const auto* col : columns)
            std::cout << std::setw(colWidth) << col->getValue(i);
        std::cout << "\n";
    }
    std::cout << std::defaultfloat;
}
