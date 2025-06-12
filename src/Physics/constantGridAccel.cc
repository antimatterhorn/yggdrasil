// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include <iostream>

template <int dim>
class ConstantGridAccel : public Physics<dim> {
protected:
    Lin::Vector<dim> gravityVector;
    double dtmin;

public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    ConstantGridAccel(NodeList* nodeList, PhysicalConstants& constants, Vector& gravityVector) :
        Physics<dim>(nodeList, constants),
        gravityVector(gravityVector) {

        int numZones = nodeList->size();
        this->template EnrollFields<Vector>({"acceleration", "velocity", "position"});

        for (int i = 0; i < numZones; ++i)
            nodeList->getField<Vector>("acceleration")->setValue(i, gravityVector);

        this->template EnrollStateFields<Vector>({"velocity"});
    }

    virtual void EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        int numZones = nodeList->size();

        VectorField* acceleration = nodeList->getField<Vector>("acceleration");
        VectorField* velocity     = initialState->template getField<Vector>("velocity");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        dvdt->copyValues(acceleration);  // dv/dt = a

        this->lastDt = dt;
    }

    virtual double EstimateTimestep() const override {
        return 1e30; // this physics package does not support setting the timestep for now
    }

    virtual std::string name() const override { return "constantGridAccel"; }
    virtual std::string description() const override {
        return "Constant acceleration in the grid"; }
};
