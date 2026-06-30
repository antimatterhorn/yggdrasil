// Copyright (C) 2025  Cody Raskin

#include "physics.hh"

// Applies a constant body acceleration to every particle in a DEM NodeList.
//
// Designed as a companion package to DEM<dim>.  Because the integrator drives
// each package independently and writes back to the shared NodeList in sequence,
// two packages that both update position would double-advance it.  This class
// avoids that by enrolling ONLY velocity as a state field: the default
// FinalizeStep only copies velocity back to the NodeList; DEM owns position.
//
// Place this package BEFORE DEM in the packages list so that the body-force
// velocity kick is applied before DEM reads the state:
//
//   packages = [DEMConstantForce2d(...), DEM2d(...)]

template <int dim>
class DEMConstantForce : public Physics<dim> {
protected:
    Lin::Vector<dim> bodyForce;   // constant acceleration [length / time²]

public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    DEMConstantForce(NodeList* nodeList, PhysicalConstants& constants, Vector bodyForce) :
        Physics<dim>(nodeList, constants),
        bodyForce(bodyForce) {
        // Enroll only velocity — position is owned by DEM
        this->template EnrollStateFields<Vector>({"velocity"});
    }

    ~DEMConstantForce() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        int n = this->nodeList->size();
        VectorField* dvdt = deriv.template getField<Vector>("velocity");
        for (int i = 0; i < n; ++i)
            dvdt->setValue(i, bodyForce);
        this->lastDt = dt;
    }

    // EstimateTimestep returns 0 so this package never constrains the global dt
    virtual double EstimateTimestep() const override { return 0.0; }

    virtual std::string name() const override { return "DEMConstantForce"; }
    virtual std::string description() const override {
        return "Constant body acceleration companion package for DEM";
    }
};
