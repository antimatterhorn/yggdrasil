// Copyright (C) 2025  Cody Raskin

#include "kinematics.hh"
#include <iostream>

template <int dim>
class NBodyGravity : public Kinematics<dim> {
protected:
    double dtmin;
    double plummerLength;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    NBodyGravity(NodeList* nodeList, PhysicalConstants& constants, double plummerLength) :
        Kinematics<dim>(nodeList,constants),
        plummerLength(plummerLength) {}

    ~NBodyGravity() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        ScalarField* mass           = nodeList->getField<double>("mass");
        VectorField* position       = initialState->template getField<Vector>("position");
        VectorField* acceleration   = nodeList->getField<Vector>("acceleration");
        VectorField* velocity       = initialState->template getField<Vector>("velocity");

        VectorField* dxdt           = deriv.template getField<Vector>("position");
        VectorField* dvdt           = deriv.template getField<Vector>("velocity");

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i=0; i<numNodes ; ++i) {
            Vector a = Vector::zero();
            Vector rij = Vector();
            for (int j=0; j<numNodes; ++j)
            {
                if (j!=i) {
                    rij = position->getValue(j) - position->getValue(i);
                    double mj = mass->getValue(j);
                    a += constants.G()*mj/(rij.mag2()+plummerLength)*rij.normal();
                }
            }
            Vector v = velocity->getValue(i);
            acceleration->setValue(i,a);
            double amag = a.mag2();
            double vmag = v.mag2();
            local_dtmin = std::min(local_dtmin,vmag/amag);
            dxdt->setValue(i,v+dt*a);
            dvdt->setValue(i,a);
        }
        dtmin  = local_dtmin;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 1e-2; // Adjust as needed
        double timestep = timestepCoefficient * sqrt(dtmin);

        return timestep;
    }

    virtual std::string name() const override { return "nBodyGravity"; }
    virtual std::string description() const override {
        return "N-body gravity physics package for particles"; }
};
