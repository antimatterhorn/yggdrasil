// Copyright (C) 2026  Cody Raskin

#pragma once

#include <cmath>
#include "../Math/vectorMath.hh"
#include "../Math/tensorMath.hh"

// Derivative tensors of the Newtonian kernel f(R) = 1/|R|, dimension-generic
namespace FMM {

// D1(R) = grad(1/|R|) = -R / |R|^3
template <int dim>
Lin::Vector<dim> D1(const Lin::Vector<dim>& R) {
    double r = R.magnitude();
    return R * (-1.0 / (r*r*r));
}

// D2(R) = grad grad (1/|R|) = (3 R⊗R - |R|^2 I) / |R|^5
template <int dim>
Lin::Tensor<dim> D2(const Lin::Vector<dim>& R) {
    double r2 = R.mag2();
    double r5 = r2*r2*std::sqrt(r2);
    Lin::Tensor<dim> RR = Lin::Tensor<dim>::outer(R, R);
    return (RR * 3.0 - Lin::Tensor<dim>::one() * r2) * (1.0 / r5);
}

// (M : D3)_c = sum_ab M_ab D3_abc, for symmetric rank-2 M. Returns a vector.
template <int dim>
Lin::Vector<dim> quadrupoleGradTerm(const Lin::Tensor<dim>& M, const Lin::Vector<dim>& R) {
    double r2 = R.mag2();
    double r5 = r2*r2*std::sqrt(r2);
    double r7 = r5*r2;
    double t = M.trace();
    Lin::Vector<dim> u = M.vecMul(R);
    double q = R.dot(u);
    return (R*(3.0*t) + u*6.0) * (1.0/r5) - R * (15.0*q/r7);
}

// (M : D4)_cd = sum_ab M_ab D4_abcd, for symmetric rank-2 M. Returns a rank-2 tensor.
template <int dim>
Lin::Tensor<dim> quadrupoleHessTerm(const Lin::Tensor<dim>& M, const Lin::Vector<dim>& R) {
    double r2 = R.mag2();
    double r5 = r2*r2*std::sqrt(r2);
    double r7 = r5*r2;
    double r9 = r7*r2;
    double t = M.trace();
    Lin::Vector<dim> u = M.vecMul(R);
    double q = R.dot(u);

    Lin::Tensor<dim> I  = Lin::Tensor<dim>::one();
    Lin::Tensor<dim> RR = Lin::Tensor<dim>::outer(R, R);
    Lin::Tensor<dim> uR = Lin::Tensor<dim>::outer(u, R);
    Lin::Tensor<dim> Ru = Lin::Tensor<dim>::outer(R, u);

    Lin::Tensor<dim> term1 = (I*(3.0*t) + M*6.0) * (1.0/r5);
    Lin::Tensor<dim> term2 = (RR*t + I*q) * (-15.0/r7);
    Lin::Tensor<dim> term3 = (uR + Ru) * (-30.0/r7);
    Lin::Tensor<dim> term4 = RR * (105.0*q/r9);
    return term1 + term2 + term3 + term4;
}

} // namespace FMM
