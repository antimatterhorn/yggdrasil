// Copyright (C) 2025  Cody Raskin

#include "integrator.hh"

template <int dim>
class RungeKutta2Integrator : public Integrator<dim> {
protected:

public:
    RungeKutta2Integrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose = false) : 
        Integrator<dim>(packages,dtmin,verbose) {}

    ~RungeKutta2Integrator() {}

    virtual State<dim> 
    Integrate(Physics<dim>* physics) override {
        double dt = this->dt;
        double time = this->time;

        const State<dim>* state  = physics->getState();
        State<dim> interim = state->deepCopy();
        State<dim> k1(state->size());
        State<dim> k2(state->size());

        k1.ghost(state);
        k2.ghost(state);

        physics->EvaluateDerivatives(state,k1,time,0);

        interim +=k1*dt; 

        physics->EvaluateDerivatives(&interim,k2,time,dt);

        State<dim> newState = state->deepCopy();

        k1 += k2;
        k1 *= 0.5*dt;
        newState += k1;

        return newState;
    }
};
