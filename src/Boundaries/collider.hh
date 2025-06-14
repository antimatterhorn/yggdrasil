#pragma once

#include <vector>
#include "../Physics/physics.hh"
#include "boundary.hh"

// Base class for Grid Boundaries
template <int dim>
class Collider : public Boundary<dim> {
protected:

public:
    Collider(){}
    
    virtual ~Collider() {}

    virtual bool Inside(Lin::Vector<dim>& otherPosition, double otherRadius) { return false; }

};