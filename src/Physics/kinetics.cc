// Copyright (C) 2025  Cody Raskin

#include "kinematics.hh"
#include <iostream>

template <int dim>
class Kinetics : public Kinematics<dim> {
protected:
    double dtmin;
    double timeVisited;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Kinetics(NodeList* nodeList, PhysicalConstants& constants) :
        Kinematics<dim>(nodeList,constants) {

        this->template EnrollFields<double>({"mass", "radius"});

        // kinetics will operate solely on the spatial derivative. velocities merely change direction
        // i.e. there is no acceleration term for this physics package so no dvdt
    }

    ~Kinetics() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {
        if (time != timeVisited) {
            PhysicalConstants constants = this->constants;
            int numNodes = this->nodeList->size();

            ScalarField* mass           = this->nodeList->template getField<double>("mass");
            ScalarField* radius         = this->nodeList->template getField<double>("radius");
            VectorField* position       = initialState->template getField<Vector>("position");
            VectorField* velocity       = this->nodeList->template getField<Vector>("velocity");

            //VectorField* dxdt           = deriv.template getField<Vector>("position");
            double local_dtmin = 1e30;

            #pragma omp parallel for reduction(min:local_dtmin)
            for (int i = 0; i < numNodes - 1; ++i) {
                dtmin = 1e30;
                for (int j = i + 1; j < numNodes; ++j) {
                    Vector ri = position->getValue(i);
                    Vector rj = position->getValue(j);
                    double si = radius->getValue(i);
                    double sj = radius->getValue(j);
                    Vector vi = velocity->getValue(i);
                    Vector vj = velocity->getValue(j);

                    Vector vij = vi-vj;
                    Vector rij = ri-rj;

                    if (rij.magnitude() < (si+sj) && vij*rij<0) { //collision
                        double mi = mass->getValue(i);
                        double mj = mass->getValue(j);

                        Vector vip = vi - 2*mj/(mi+mj)*(vij*rij/rij.mag2())*rij;
                        Vector vjp = vj + 2*mi/(mi+mj)*(vij*rij/rij.mag2())*rij;

                        velocity->setValue(i,vip);
                        velocity->setValue(j,vjp);
                        // dxdt->setValue(i,vip);
                        // dxdt->setValue(j,vjp);

                        continue;  // skip this i (1 collision only please!)
                    }
                }
                Vector vi = velocity->getValue(i);
                double si = radius->getValue(i);
                local_dtmin = std::min(local_dtmin,0.25*si/vi.magnitude());
            }
            dtmin = local_dtmin;
        }
        
        
        timeVisited = time;
        this->lastDt = dt;
    }

    virtual void
    FinalizeStep(const State<dim>* finalState) override {
        // don't do anything with the state since kinetics modifies nodelist directly
    };

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 0.25; // Adjust as needed
        double timestep = timestepCoefficient * sqrt(dtmin);

        return timestep;
    }

    virtual std::string name() const override { return "kinetics"; }
    virtual std::string description() const override {
        return "The kinetics physics package for particles"; }
};
