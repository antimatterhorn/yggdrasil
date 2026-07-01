// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include <iostream>

template <int dim>
class ConstantGravity : public Physics<dim> {
protected:
    Lin::Vector<dim> gravityVector;
    double dtmin;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    ConstantGravity(NodeList* nodeList, PhysicalConstants& constants, Vector& gravityVector) :
        Physics<dim>(nodeList,constants),
        gravityVector(gravityVector) {

        this->template EnrollFields<Vector>({"acceleration"});
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);

        int numNodes = nodeList->size();
        for (int i=0; i<numNodes; ++i)
            nodeList->getField<Vector>("acceleration")->setValue(i, gravityVector);
    }

    ~ConstantGravity() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        int numNodes = nodeList->size();

        VectorField* acceleration = nodeList->getField<Vector>("acceleration");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        dvdt->copyValues(acceleration);

        double local_dtmin = 1e30;
        VectorField* velocity = initialState->template getField<Vector>("velocity");

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i=0; i<numNodes; ++i) {
            double amag = acceleration->getValue(i).mag2();
            double vmag = velocity->getValue(i).mag2();
            if (amag > 0.0)
                local_dtmin = std::min(local_dtmin, vmag / amag);
        }
        dtmin = local_dtmin;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 1e-4;
        return timestepCoefficient * sqrt(dtmin);
    }

    virtual std::string name() const override { return "constantGravity"; }
    virtual std::string description() const override {
        return "Constant acceleration"; }
};
