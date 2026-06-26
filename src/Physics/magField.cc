// Copyright (C) 2025  Cody Raskin
// this method is intended for probe particles in a mag field - it will not work for grid-based physics

#include "kinematics.hh"
#include <iostream>

template <int dim>
class MagField : public Kinematics<dim> {
private:
    Lin::Vector<dim> B;
    double dtmin;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    MagField(NodeList* nodeList, PhysicalConstants& constants, Vector& B) :
        Kinematics<dim>(nodeList,constants),
        B(B) {
            
        this->template EnrollFields<double>({"charge"});
    }

    ~MagField() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        NodeList* nodeList = this->nodeList;
        PhysicalConstants constants = this->constants;
        int numNodes = nodeList->size();

        VectorField* velocity     = initialState->template getField<Vector>("velocity");
        ScalarField* charge       = nodeList->getField<double>("charge");
        ScalarField* mass         = nodeList->getField<double>("mass");

        VectorField* dxdt         = deriv.template getField<Vector>("position");
        VectorField* dvdt         = deriv.template getField<Vector>("velocity");

        dxdt->copyValues(velocity);

        dtmin = 1e30;

        #pragma omp parallel for reduction(min:dtmin)
        for (int i = 0; i < numNodes; ++i) {
            double q = charge->getValue(i);
            double m = mass->getValue(i);
            Vector v = velocity->getValue(i);

            Vector a = (q / m) * v.crossProduct(B);
            dvdt->setValue(i, a);

            double amag = a.mag2();
            double vmag = v.mag2();
            if (amag > 1e-10)
                dtmin = std::min(dtmin, vmag / amag);
        }

        this->lastDt = dt;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 5e-1; // Adjust as needed
        double timestep = timestepCoefficient * dtmin;

        return timestep;
    }

    void
    SetB(Vector& mF) { B = mF; }

    virtual std::string name() const override { return "magField"; }
    virtual std::string description() const override {
        return "Magnetic Field for Kinematics"; }
};
