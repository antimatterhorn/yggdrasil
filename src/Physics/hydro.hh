// Copyright (C) 2025  Cody Raskin

#ifndef HYDRO_HH
#define HYDRO_HH

#include "physics.hh"
#include "../EOS/equationOfState.hh"

template <int dim>
class Hydro : public Physics<dim> {
protected:
    EquationOfState* eos;
public:

    Hydro(NodeList* nodeList, PhysicalConstants& constants, EquationOfState* eos) : 
        Physics<dim>(nodeList,constants), eos(eos) {
        this->template EnrollFields<double>({"pressure", "density", "specificInternalEnergy", "soundSpeed"});
    }

    virtual ~Hydro() {}

    virtual std::string name() const override { return "hydro"; }
    virtual std::string description() const override {
        return "Some kind of hydro"; }
};

#endif