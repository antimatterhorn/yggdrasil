// Copyright (C) 2025  Cody Raskin
// this method is intended for probe particles in a mag field - it will not work for grid-based physics

#include "kinematics.hh"
#include <iostream>

template <int dim>
class MagField : public Kinematics<dim> {
private:
    Lin::Vector<dim> magF;
    double dtmin;
public:
    using Vector = Lin::Vector<dim>;
    using VectorField = Field<Vector>;
    using ScalarField = Field<double>;

    MagField(NodeList* nodeList, PhysicalConstants& constants, Vector& magF) :
        Kinematics<dim>(nodeList,constants),
        magF(magF) {
            
        this->template EnrollFields<double>({"charge"});
    }

    ~MagField() {}

    virtual void
    EvaluateDerivatives(const State<dim>* initialState, State<dim>& deriv, const double time, const double dt) override {

    }

    virtual double
    EstimateTimestep() const override {
        double timestepCoefficient = 1e-4; // Adjust as needed
        double timestep = timestepCoefficient * sqrt(dtmin);

        return timestep;
    }

    void
    SetMagField(Vector& mF) { magF = mf; }

    virtual std::string name() const override { return "magField"; }
    virtual std::string description() const override {
        return "Magnetic Field for Kinematics"; }
};
