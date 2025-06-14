#pragma once

#include <vector>
#include "../Physics/physics.hh"

template <int dim>
class Physics;

// Base class for Grid Boundaries
template <int dim>
class Boundary {
protected:

public:
    Boundary() {}
    
    virtual ~Boundary() {}

    virtual void ZeroTimeInitialize() {}

    virtual void
    ApplyBoundaries(State<dim>* state, NodeList* nodeList) {}

};