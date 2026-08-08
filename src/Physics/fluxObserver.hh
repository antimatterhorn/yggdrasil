// Copyright (C) 2026  Cody Raskin

#pragma once
#include "../Math/vectorMath.hh"

// Optional tap on a finite-volume solver's face fluxes, so an AMR flux register
// can accumulate them without Physics depending on the AMR module. Flux
// components are passed directly rather than as an HLLFlux to keep this
// interface free of the solver's own types.
template <int dim>
class FluxObserver {
public:
    virtual ~FluxObserver() {}

    // `plusSide` distinguishes the face on the +axis side of `cell` from the one on its -axis side.
    virtual void recordFlux(int cell, int axis, bool plusSide,
                            double massFlux,
                            const Lin::Vector<dim>& momentumFlux,
                            double energyFlux,
                            double area) = 0;
};
