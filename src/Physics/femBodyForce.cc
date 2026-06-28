// Copyright (C) 2025  Cody Raskin

#include "physics.hh"

// Abstract base for external body force packages used with FEM.
//
// Inherits Physics directly (not Kinematics) so it owns no position/velocity
// state fields and does not advance kinematics. Subclasses implement
// EvaluateDerivatives, which writes per-node accelerations into the shared
// "acceleration" NodeList field. FEMLinearElasticity reads that field and
// adds it to the elastic acceleration each step.
//
// Ordering in the integrator: add BodyForce packages BEFORE FEMLinearElasticity
// so the field is populated before FEM reads it.

template <int dim>
class FEMBodyForce : public Physics<dim> {
public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;

    FEMBodyForce(NodeList* nodeList, PhysicalConstants& constants)
        : Physics<dim>(nodeList, constants) {
        this->template EnrollFields<Vector>({"acceleration"});
    }

    // Subclasses write body force accelerations into the "acceleration" field.
    // time is the start-of-step time; dt is the current timestep (0 for the
    // first RK stage).  No derivatives are written — the deriv state is empty.
    void EvaluateDerivatives(const State<dim>* state, State<dim>& deriv,
                             const double time, const double dt) override = 0;

    double EstimateTimestep() const override { return 1e30; }

    std::string name()        const override { return "FEMBodyForce"; }
    std::string description() const override { return "Abstract body force package for FEM"; }
};

// Uniform constant body force (e.g. gravity) applied to all nodes.
template <int dim>
class ConstantBodyForce : public FEMBodyForce<dim> {
    Lin::Vector<dim> g;
public:
    using Vector      = Lin::Vector<dim>;
    using VectorField = Field<Vector>;

    ConstantBodyForce(NodeList* nodeList, PhysicalConstants& constants, Vector g)
        : FEMBodyForce<dim>(nodeList, constants), g(g) {}

    void EvaluateDerivatives(const State<dim>* state, State<dim>& deriv,
                             const double time, const double dt) override {
        VectorField* accelField = this->nodeList->template getField<Vector>("acceleration");
        int numNodes = this->nodeList->size();
        for (int i = 0; i < numNodes; ++i)
            accelField->setValue(i, g);
    }

    std::string name()        const override { return "ConstantBodyForce"; }
    std::string description() const override {
        return "Constant uniform body force applied to all nodes";
    }
};
