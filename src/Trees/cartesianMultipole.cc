// Copyright (C) 2026  Cody Raskin

#ifndef CARTESIANMULTIPOLE_CC
#define CARTESIANMULTIPOLE_CC

#include "cartesianMultipole.hh"

template <int dim>
CartesianMultipole<dim> CartesianMultipole<dim>::fromSources(const std::vector<Vector>& positions,
                                                               const std::vector<double>& masses) {
    CartesianMultipole<dim> result;
    Vector centroid = Vector::zero();
    for (size_t i = 0; i < positions.size(); ++i) {
        result.mass += masses[i];
        centroid += positions[i] * masses[i];
    }
    centroid = centroid * (1.0 / result.mass);
    result.center = centroid;

    for (size_t i = 0; i < positions.size(); ++i) {
        Vector rho = positions[i] - centroid;
        result.quadrupole = result.quadrupole + Tensor::outer(rho, rho) * masses[i];
    }
    return result;
}

template <int dim>
CartesianMultipole<dim> CartesianMultipole<dim>::combine(const std::vector<const CartesianMultipole<dim>*>& children) {
    CartesianMultipole<dim> result;
    Vector centroid = Vector::zero();
    for (const auto* child : children) {
        result.mass += child->mass;
        centroid += child->center * child->mass;
    }
    centroid = centroid * (1.0 / result.mass);
    result.center = centroid;

    for (const auto* child : children) {
        Vector d = child->center - centroid;
        result.quadrupole = result.quadrupole + child->quadrupole + Tensor::outer(d, d) * child->mass;
    }
    return result;
}

template <int dim>
LocalExpansion<dim> LocalExpansion<dim>::fromMultipole(const Vector& targetCenter,
                                                         const CartesianMultipole<dim>& source,
                                                         double G) {
    LocalExpansion<dim> result;
    result.center = targetCenter;

    Vector R0 = targetCenter - source.center;
    double r0 = R0.magnitude();
    const Tensor& M2 = source.quadrupole;
    double M0 = source.mass;

    double t0 = M2.trace();
    Vector u0 = M2.vecMul(R0);
    double q0 = R0.dot(u0);
    double quadScalar = (3.0*q0 - t0*r0*r0) / (r0*r0*r0*r0*r0);

    result.phi0 = -G * (M0/r0 + 0.5*quadScalar);
    result.g    = (FMM::D1<dim>(R0) * M0 + FMM::quadrupoleGradTerm<dim>(M2, R0) * 0.5) * (-G);
    result.H    = (FMM::D2<dim>(R0) * M0 + FMM::quadrupoleHessTerm<dim>(M2, R0) * 0.5) * (-G);

    return result;
}

template <int dim>
LocalExpansion<dim>& LocalExpansion<dim>::operator+=(const LocalExpansion<dim>& other) {
    phi0 += other.phi0;
    g    += other.g;
    H    = H + other.H;
    return *this;
}

template <int dim>
LocalExpansion<dim> LocalExpansion<dim>::shift(const Vector& newCenter) const {
    LocalExpansion<dim> result;
    result.center = newCenter;

    Vector d = center - newCenter;
    result.phi0 = phi0 - g.dot(d) + 0.5 * d.dot(H.vecMul(d));
    result.g    = g - H.vecMul(d);
    result.H    = H;
    return result;
}

template <int dim>
double LocalExpansion<dim>::evaluatePotential(const Vector& point) const {
    Vector u = point - center;
    return phi0 + g.dot(u) + 0.5 * u.dot(H.vecMul(u));
}

template <int dim>
typename LocalExpansion<dim>::Vector LocalExpansion<dim>::evaluateForce(const Vector& point) const {
    Vector u = point - center;
    return -g - H.vecMul(u);
}

#endif
