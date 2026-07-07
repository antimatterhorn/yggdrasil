// Copyright (C) 2026  Cody Raskin

#pragma once

// Non-templated interface exposing the parts of Integrator<dim> that don't
// actually depend on dim: the simulation clock. Code that only needs to
// read or restore cycle/time/dt -- e.g. IO::RestartWriter/RestartReader --
// can hold an IntegratorBase& instead of an Integrator<dim>&, so it doesn't
// need a dimension-specific variant of its own.
class IntegratorBase {
public:
    virtual ~IntegratorBase() {}

    virtual double const Time() = 0;
    virtual unsigned int Cycle() = 0;
    virtual double const Dt() = 0;
    virtual void restoreState(unsigned int cycle, double time, double dt) = 0;
};
