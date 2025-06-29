// Copyright (C) 2025  Cody Raskin

#pragma once

#include "integrator.hh"

template <int dim>
Integrator<dim>::Integrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose)
    : packages(packages), dt(dtmin), dtmin(dtmin), cycle(0), time(0.0), verbose(verbose) {}

template <int dim>
Integrator<dim>::~Integrator() {}

template <int dim>
void Integrator<dim>::Step() {  
    if (cycle == 0) {
        for (Physics<dim>* physics : packages)
            physics->ZeroTimeInitialize();
    }

    for (Physics<dim>* physics : packages)
    {
        physics->UpdateState();
        physics->PreStepInitialize();
        
        State<dim> finalState = Integrate(physics);

        physics->ApplyBoundaries(&finalState);
        physics->FinalizeStep(&finalState);
    }

    time += dt;
    cycle += 1;

    VoteDt();
}

template <int dim>
State<dim> 
Integrator<dim>::Integrate(Physics<dim>* physics) {
    const State<dim>* state = physics->getState();
    State<dim> derivatives(state->size());
    derivatives.ghost(state);

    physics->EvaluateDerivatives(state, derivatives, time, 0);

    derivatives *= dt;
    State<dim> newState(state->size());
    newState = state->deepCopy();
    newState += derivatives;
    return newState;
}

template <int dim>
void Integrator<dim>::VoteDt() {
        double smallestDt = 1e30;
        for (Physics<dim>* physics : packages) {
            double newdt = physics->EstimateTimestep();
            if (newdt < smallestDt) {
                smallestDt = newdt;
                if (verbose)
                    std::cout << physics->name() << " requested timestep of " << newdt << "\n";
            }
        }

        dt = (dt < smallestDt ?  dt + 0.2 * (smallestDt - dt) : smallestDt);

        this->dt = std::max(dt, this->dtmin) * this->dtMultiplier;
}

template <int dim>
double const Integrator<dim>::Time() { return time; }

template <int dim>
unsigned int Integrator<dim>::Cycle() { return cycle; }

template <int dim>
double const Integrator<dim>::Dt() { return dt; }
