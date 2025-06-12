// Copyright (C) 2025  Cody Raskin

#include "kinematics.hh"
#include <iostream>

template <int dim>
class ConstantGravity : public Kinematics<dim> {
protected:
    Lin::Vector<dim> gravityVector;
    double dtmin;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    ConstantGravity(NodeList* nodeList, PhysicalConstants& constants, Vector& gravityVector) :
        Kinematics<dim>(nodeList,constants),
        gravityVector(gravityVector) {

        int numNodes = nodeList->size();
        for (int i=0; i<numNodes; ++i)
            nodeList->getField<Vector>("acceleration")->setValue(i,gravityVector);
    }

    ~ConstantGravity() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        VectorField* position       = initialState->template getField<Vector>("position");
        VectorField* acceleration   = nodeList->getField<Vector>("acceleration");
        // since the acceleration is literally constant, we keep it on the nodeList
        VectorField* velocity       = initialState->template getField<Vector>("velocity");

        VectorField* dxdt           = deriv.template getField<Vector>("position");
        VectorField* dvdt           = deriv.template getField<Vector>("velocity");

        dxdt->copyValues(velocity);
        dxdt->operator+(*acceleration*dt);
        dvdt->copyValues(acceleration);

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i=0; i<numNodes ; ++i) {
            double amag = acceleration->getValue(i).mag2();
            double vmag = velocity->getValue(i).mag2();
            local_dtmin = std::min(local_dtmin,vmag/amag);
        }
        dtmin = local_dtmin;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 1e-4; // Adjust as needed
        double timestep = timestepCoefficient * sqrt(dtmin);

        return timestep;
    }

    virtual std::string name() const override { return "constantGravity"; }
    virtual std::string description() const override {
        return "Constant acceleration"; }
};
