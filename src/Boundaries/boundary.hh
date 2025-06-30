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
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;
    using Complex     = std::complex<double>;
    using ComplexField= Field<Complex>;
    
    Boundary() {}
    
    virtual ~Boundary() {}

    virtual void ZeroTimeInitialize() {}

    virtual void
    ApplyBoundaries(State<dim>* state, NodeList* nodeList) {
        for (int i = 0; i < state->count(); ++i) {
            FieldBase* field = state->getFieldByIndex(i);
            if (typeid(*field) == typeid(ScalarField)) {
                ApplyThis(static_cast<ScalarField*>(field));
            } else if (typeid(*field) == typeid(VectorField)) {
                ApplyThis(static_cast<VectorField*>(field));
            } else if (typeid(*field) == typeid(ComplexField)) {
                ApplyThis(static_cast<ComplexField*>(field));
            }
        }
    }

    virtual void 
    ApplyThis(ScalarField* field) {};
    virtual void
    ApplyThis(VectorField* field) {};
    virtual void
    ApplyThis(ComplexField* field) {};

};