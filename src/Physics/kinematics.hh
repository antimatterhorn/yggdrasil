// Copyright (C) 2025  Cody Raskin

#pragma once

#include "physics.hh"
#include <iostream>

template <int dim>
class Kinematics : public Physics<dim> {
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    Kinematics(NodeList* nodeList, PhysicalConstants& constants) :
        Physics<dim>(nodeList,constants) {
            
        this->template EnrollFields<double>({"mass"});
        this->template EnrollFields<Vector>({"acceleration", "velocity", "position"});
        this->template EnrollStateFields<Vector>({"velocity", "position"}, FieldPolicy::INTEGRATE);
    }

    ~Kinematics() {}

    // Transport step: advance positions by velocity. Force packages add to dvdt via ACCUMULATE.
    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv,
                        const double time, const double dt) override {
        VectorField* velocity = initialState->template getField<Vector>("velocity");
        VectorField* dxdt     = deriv.template getField<Vector>("position");
        dxdt->copyValues(velocity);
    }

    virtual std::string name() const override { return "kinematics"; }
    virtual std::string description() const override {
        return "Kinematics physics package for particles"; }
};
