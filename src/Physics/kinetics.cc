// Copyright (C) 2026  Cody Raskin

#include "kinematics.hh"
#include <iostream>

template <int dim>
class Kinetics : public Kinematics<dim> {
protected:
    double dtmin = 1e30;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Kinetics(NodeList* nodeList, PhysicalConstants& constants) :
        Kinematics<dim>(nodeList,constants) {

        this->template EnrollFields<double>({"mass", "radius"});

        // Kinematics base provides EvaluateDerivatives (dxdt = v).
        // Force packages ACCUMULATE dvdt on top of that.
        // Collision response happens post-integration in FinalChecks.
    }

    ~Kinetics() {}

    // Collision detection and response run after FinalizeStep has written the
    // new positions/velocities to the NodeList, so we work on NodeList directly.
    virtual void
    FinalChecks() override {
        NodeList* nodeList = this->nodeList;
        int numNodes = nodeList->size();

        ScalarField* mass     = nodeList->template getField<double>("mass");
        ScalarField* radius   = nodeList->template getField<double>("radius");
        VectorField* position = nodeList->template getField<Vector>("position");
        VectorField* velocity = nodeList->template getField<Vector>("velocity");

        double local_dtmin = 1e30;

        for (int i = 0; i < numNodes - 1; ++i) {
            for (int j = i + 1; j < numNodes; ++j) {
                Vector ri = position->getValue(i);
                Vector rj = position->getValue(j);
                double si = radius->getValue(i);
                double sj = radius->getValue(j);
                Vector vi = velocity->getValue(i);
                Vector vj = velocity->getValue(j);

                Vector vij = vi - vj;
                Vector rij = ri - rj;

                if (rij.magnitude() < (si + sj) && vij * rij < 0) {
                    double mi = mass->getValue(i);
                    double mj = mass->getValue(j);

                    Vector vip = vi - 2*mj/(mi+mj) * (vij*rij/rij.mag2()) * rij;
                    Vector vjp = vj + 2*mi/(mi+mj) * (vij*rij/rij.mag2()) * rij;

                    velocity->setValue(i, vip);
                    velocity->setValue(j, vjp);
                }
            }

            Vector vi = velocity->getValue(i);
            double si = radius->getValue(i);
            if (vi.magnitude() > 0.0)
                local_dtmin = std::min(local_dtmin, si/vi.magnitude());
        }

        dtmin = local_dtmin;
    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 0.2;
        return timestepCoefficient * dtmin;
    }

    virtual std::string name() const override { return "kinetics"; }
    virtual std::string description() const override {
        return "The kinetics physics package for particles"; }
};
