// Copyright (C) 2025  Cody Raskin

#ifndef GRIDBoundary_HH
#define GRIDBoundary_HH

#include <vector>
#include "../Mesh/grid.hh"
#include "../Physics/physics.hh"
#include "boundary.hh"

// Base class for Grid Boundary
template <int dim>
class GridBoundary : public Boundary<dim> {
protected:
    Mesh::Grid<dim>* grid;
    std::vector<bool> activeFaces;

public:
    GridBoundary(Mesh::Grid<dim>* grid) :
        grid(grid), activeFaces(2 * dim, true) {}

    virtual ~GridBoundary() {}

    virtual void setFaces(const std::vector<std::string>& faces) {
        static const std::array<std::string, 6> names =
            {"left", "right", "bottom", "top", "front", "back"};
        std::fill(activeFaces.begin(), activeFaces.end(), false);
        for (const auto& f : faces)
            for (int i = 0; i < 2 * dim; ++i)
                if (names[i] == f) { activeFaces[i] = true; break; }
    }
};

#endif // GRIDBoundary_HH