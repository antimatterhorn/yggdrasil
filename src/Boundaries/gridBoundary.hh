// Copyright (C) 2026  Cody Raskin

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

    bool refreshPerStage() const override { return true; }

    // Face names in [-x, +x, -y, +y, -z, +z] order, so names[i] labels the same
    // face that Grid::lowMost(i/2) / highMost(i/2) returns and that every
    // concrete subclass stores at index i.
    virtual void setFaces(const std::vector<std::string>& faces) {
        static const std::array<std::string, 6> names =
            {"left", "right", "bottom", "top", "front", "back"};
        std::fill(activeFaces.begin(), activeFaces.end(), false);
        for (const auto& f : faces)
            for (int i = 0; i < 2 * dim; ++i)
                if (names[i] == f) { activeFaces[i] = true; break; }
    }

    // Whether face `i` (in [-x,+x,-y,+y,-z,+z] order) is active.
    bool isFaceActive(int i) const { return activeFaces[i]; }
};

#endif // GRIDBoundary_HH