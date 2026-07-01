// Copyright (C) 2025  Cody Raskin

#include "physics.hh"
#include <iostream>

template <int dim>
class NBodyGravity : public Physics<dim> {
protected:
    double dtmin;
    double plummerLength;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    NBodyGravity(NodeList* nodeList, PhysicalConstants& constants, double plummerLength) :
        Physics<dim>(nodeList,constants),
        plummerLength(plummerLength) {

        this->template EnrollFields<Vector>({"acceleration"});
        // ACCUMULATE: contributes dvdt but does not own or finalize velocity.
        // position is read from the INTEGRATE package's sub-stage state via initialState.
        this->template EnrollStateFields<Vector>({"velocity"}, FieldPolicy::ACCUMULATE);
    }

    ~NBodyGravity() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        ScalarField* mass         = nodeList->getField<double>("mass");
        VectorField* position     = initialState->template getField<Vector>("position");
        VectorField* velocity     = initialState->template getField<Vector>("velocity");
        VectorField* acceleration = nodeList->getField<Vector>("acceleration");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        double local_dtmin = 1e30;

        #pragma omp parallel for reduction(min:local_dtmin)
        for (int i=0; i<numNodes; ++i) {
            Vector a = Vector::zero();
            for (int j=0; j<numNodes; ++j) {
                if (j != i) {
                    Vector rij = position->getValue(j) - position->getValue(i);
                    double mj = mass->getValue(j);
                    a += constants.G() * mj / (rij.mag2() + plummerLength) * rij.normal();
                }
            }
            Vector v = velocity->getValue(i);
            acceleration->setValue(i, a);
            dvdt->setValue(i, a);

            double amag = a.mag2();
            double vmag = v.mag2();
            if (amag > 0.0)
                local_dtmin = std::min(local_dtmin, vmag / amag);
        }
        dtmin = local_dtmin;
        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 1e-2;
        return timestepCoefficient * sqrt(dtmin);
    }

    virtual std::string name() const override { return "nBodyGravity"; }
    virtual std::string description() const override {
        return "N-body gravity physics package for particles"; }
};
