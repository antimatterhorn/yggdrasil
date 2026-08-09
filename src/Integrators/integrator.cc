// Copyright (C) 2026  Cody Raskin

#pragma once

#include "integrator.hh"

template <int dim>
Integrator<dim>::Integrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose)
    : packages(packages), dt(dtmin), dtmin(dtmin), cycle(0), time(0.0), verbose(verbose) {}

template <int dim>
Integrator<dim>::~Integrator() {}

template <int dim>
void Integrator<dim>::Initialize() {
    if (initialized) return;
    for (Physics<dim>* physics : packages)
        physics->ZeroTimeInitialize();
    initialized = true;
}

template <int dim>
void Integrator<dim>::Step() {
    Initialize();

    // Snapshot all states before any FinalizeStep modifies the NodeList.
    for (Physics<dim>* physics : packages) {
        physics->UpdateState();
        physics->PreStepInitialize();
    }

    for (Physics<dim>* physics : packages) {
        if (!physics->hasIntegrateFields()) continue;

        std::vector<Physics<dim>*> accumulators = findAccumulators(physics);
        State<dim> finalState = Integrate(physics, accumulators);

        physics->ApplyBoundaries(&finalState);
        physics->FinalizeStep(&finalState);
    }

    time += dt;
    cycle += 1;

    VoteDt();
}

// Base (Euler) integration — concrete integrators override this.
template <int dim>
State<dim>
Integrator<dim>::Integrate(Physics<dim>* physics,
                            const std::vector<Physics<dim>*>& accumulators) {
    const State<dim>* state = physics->getState();
    State<dim> derivatives(state->size());
    derivatives.ghost(state);

    physics->EvaluateDerivatives(state, derivatives, time, 0);

    for (auto* acc : accumulators) {
        State<dim> accDeriv(state->size());
        accDeriv.ghost(acc->getState());
        acc->EvaluateDerivatives(state, accDeriv, time, 0);
        derivatives.accumulateFrom(accDeriv, acc->accumulateFieldNames());
    }

    State<dim> newState = state->deepCopy();
    newState.axpy(dt, derivatives);
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
std::vector<Physics<dim>*>
Integrator<dim>::findAccumulators(Physics<dim>* owner) {
    std::vector<std::string> ownedFields = owner->integrateFieldNames();
    std::vector<Physics<dim>*> result;
    for (auto* p : packages) {
        if (p == owner) continue;
        for (const std::string& name : p->accumulateFieldNames()) {
            if (std::find(ownedFields.begin(), ownedFields.end(), name) != ownedFields.end()) {
                result.push_back(p);
                break;
            }
        }
    }
    return result;
}

template <int dim>
double const Integrator<dim>::Time() { return time; }

template <int dim>
unsigned int Integrator<dim>::Cycle() { return cycle; }

template <int dim>
double const Integrator<dim>::Dt() { return dt; }

template <int dim>
void Integrator<dim>::restoreState(unsigned int cycle, double time, double dt) {
    this->cycle = cycle;
    this->time = time;
    this->dt = dt;
    // ZeroTimeInitialize() rebuilds per-process derived state (EOS lookups,
    // boundary/obstacle setup, GridHydroBase::insideIds) that a freshly
    // constructed process never ran, regardless of which cycle is restored --
    // it does not touch cycle/time/dt, so it's safe to run again here.
    this->initialized = false;
}
