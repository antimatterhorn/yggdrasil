// Copyright (C) 2026  Cody Raskin

#pragma once

#include <vector>
#include <algorithm>
#include <iostream>
#include "../Math/vectorMath.hh"
#include "../State/state.hh"
#include "../Boundaries/boundary.hh"

template <int dim>
class Physics; // forward declaration

template <int dim>
class Integrator {
protected:
    std::vector<Physics<dim>*> packages;
    unsigned int cycle;
    bool verbose;
    double time, dt, dtmin, dtMultiplier = 1;

    // Returns all packages that ACCUMULATE into any field that owner INTEGRATEs.
    std::vector<Physics<dim>*> findAccumulators(Physics<dim>* owner);

public:
    Integrator(std::vector<Physics<dim>*> packages, double dtmin, bool verbose = false);
    ~Integrator();

    virtual void Step();
    virtual State<dim> Integrate(Physics<dim>* physics,
                                  const std::vector<Physics<dim>*>& accumulators);
    virtual void VoteDt();
    virtual double const Time();
    virtual unsigned int Cycle();
    virtual double const Dt();
    inline std::vector<Physics<dim>*> getPackages() { return packages;};
};

#include "integrator.cc"
