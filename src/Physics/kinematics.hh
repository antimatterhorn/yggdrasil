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
        this->template EnrollStateFields<Vector>({"velocity", "position"});
    }

    ~Kinematics() {}

    virtual std::string name() const override { return "kinematics"; }
    virtual std::string description() const override {
        return "Kinematics physics package for particles"; }
};
