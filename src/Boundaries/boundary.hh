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
                ApplyScalar(static_cast<ScalarField*>(field));
            } else if (typeid(*field) == typeid(VectorField)) {
                ApplyVector(static_cast<VectorField*>(field));
            } else if (typeid(*field) == typeid(ComplexField)) {
                ApplyComplex(static_cast<ComplexField*>(field));
            }
        }
    }

    virtual void ApplyScalar(ScalarField* field) {};
    virtual void ApplyVector(VectorField* field) {};
    virtual void ApplyComplex(ComplexField* field) {};

};