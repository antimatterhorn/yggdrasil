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
public:
    GridBoundary(Mesh::Grid<dim>* grid) : 
        grid(grid) {}
    
    virtual ~GridBoundary() {}

};

#endif // GRIDBoundary_HH