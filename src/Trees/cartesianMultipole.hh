// Copyright (C) 2026  Cody Raskin

#pragma once

#include <vector>
#include "cartesianKernel.hh"

// Order-2 Cartesian multipole/local expansion of Phi(x) = -G*sum_i m_i/|x-x_i| 
// (NBodyGravity/TreeGravity's force law), dimension-generic
template <int dim>
class CartesianMultipole {
public:
    using Vector = Lin::Vector<dim>;
    using Tensor = Lin::Tensor<dim>;

    Vector center = Vector::zero();  // mass-weighted centroid; keeps the dipole moment identically zero
    double mass = 0.0;
    Tensor quadrupole = Tensor::zero();  // sum_i m_i * (x_i - center) (x) (x_i - center)

    // P2M: leaf multipole from a set of point masses, centered on their mass-weighted centroid.
    static CartesianMultipole<dim> fromSources(const std::vector<Vector>& positions,
                                                const std::vector<double>& masses);

    // M2M: combine multipoles already centered on their own COMs into one about their combined COM.
    static CartesianMultipole<dim> combine(const std::vector<const CartesianMultipole<dim>*>& children);
};

template <int dim>
class LocalExpansion {
public:
    using Vector = Lin::Vector<dim>;
    using Tensor = Lin::Tensor<dim>;

    Vector center = Vector::zero();
    double phi0 = 0.0;
    Vector g = Vector::zero();
    Tensor H = Tensor::zero();

    // M2L: local expansion about targetCenter of source's far-field potential; valid only if well separated.
    static LocalExpansion<dim> fromMultipole(const Vector& targetCenter,
                                              const CartesianMultipole<dim>& source,
                                              double G);

    // Sums contributions from multiple well-separated sources onto the same target center.
    LocalExpansion<dim>& operator+=(const LocalExpansion<dim>& other);

    // L2L: re-express this expansion about a new (child) center.
    LocalExpansion<dim> shift(const Vector& newCenter) const;

    // L2P
    double evaluatePotential(const Vector& point) const;
    Vector evaluateForce(const Vector& point) const;
};

#include "cartesianMultipole.cc"
