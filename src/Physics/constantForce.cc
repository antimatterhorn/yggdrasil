// Copyright (C) 2026  Cody Raskin

#include "physics.hh"
#include <iostream>

template <int dim>
class ConstantForce : public Physics<dim> {
protected:
    Lin::Vector<dim> forceVector;
    double dtmin;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    ConstantForce(NodeList* nodeList, PhysicalConstants& constants, Vector& forceVector) :
        Physics<dim>(nodeList,constants),
        forceVector(forceVector) {

        this->template EnrollFields<Vector>({"acceleration"});
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);

        int numNodes = nodeList->size();
        for (int i=0; i<numNodes; ++i)
            nodeList->getField<Vector>("acceleration")->setValue(i, forceVector);
    }

    ~ConstantForce() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        int numNodes = nodeList->size();

        VectorField* acceleration = nodeList->getField<Vector>("acceleration");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        dvdt->copyValues(acceleration);

        double amag = forceVector.mag2();
        double vmax_sq = 0.0;
        VectorField* velocity = initialState->template getField<Vector>("velocity");

        #pragma omp parallel for reduction(max:vmax_sq)
        for (int i=0; i<numNodes; ++i)
            vmax_sq = std::max(vmax_sq, velocity->getValue(i).mag2());

        dtmin = (amag > 0.0 && vmax_sq > 0.0) ? vmax_sq / amag : 1e30;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 0.5;
        return timestepCoefficient * dtmin;  // dt ~ 0.5 * v_max / |a|
    }

    virtual std::string name() const override { return "constantForce"; }
    virtual std::string description() const override {
        return "Constant acceleration"; }
};
